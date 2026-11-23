#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created 11 2026

Zarr zip store for parallel writing and reading: ZARZS.ZipStore is a drop-in for zarr.storage.ZipStore.

	store = ZARZS.ZipStore("volume.zip", mode="w")    # "r", "w" (new zip) or "a" (append)
	group = zarr.open_group(store, mode="w")
	...                                               # write from as many threads as useful
	store.close()                                     # writes the directory of the zip

The result is a standard zip with the entries stored uncompressed, as zarr's ZipStore writes them;
zarr.storage.ZipStore, Python's zipfile and unzip read it.

Why: zarr's ZipStore writes every entry through Python's zipfile under one lock, so the CRC-32 and the
write run for one chunk at a time, and it reads under the same lock. With large chunks the zip is
then much slower than a directory store. Here entries are stored uncompressed, so their size and CRC
are known before they are written: a thread computes the CRC, reserves the space of the entry under
the lock (microseconds) and writes the local header and data with os.pwrite at that offset outside
the lock. zarr calls the store on its event loop thread; reads and writes run in worker threads
(asyncio.to_thread), as for zarr's LocalStore. An entry is added to the directory of the zip after its
data is written, so readers see complete entries only; they read with os.pread. close() writes the
central directory with zipfile (zip64 where needed).

Beyond zarr's ZipStore:
	- a key written again replaces its entry: the zip has every name once (zarr's ZipStore adds a
	  duplicate entry with a warning); the bytes of the old entry stay unused in the file
	  (replaced_bytes), small for metadata such as zarr.json
	- delete() and delete_dir() remove entries (zarr's ZipStore raises NotImplementedError), so
	  arrays can be overwritten; their bytes stay unused in the file as well. As for zarr's
	  LocalStore, delete(path) of a path that is no entry deletes everything under path/ (an array
	  or group)

Appending (mode "a") is safe against a killed process: new entries are written after the end of the
file, so the zip as it was stays intact up to its original size, and that size is recorded in
<zip>.zarzs-append until close() has written the new central directory. If the process dies in
between, the zip can not be opened (the store refuses while the record exists) and
ZARZS.recover(path) truncates it to the state before the append. The old central directory stays in
the file unused (replaced_bytes counts it).

Limits: for threads of one process (copying the store to another process raises, a writer there
would truncate the zip); entries are written uncompressed (compression=zipfile.ZIP_STORED), entries of
an existing zip compressed otherwise are read through zipfile. It overrides private methods of
zarr.storage.ZipStore (zarr 3), test_ZARZS checks them.

@author: Vojtěch Kulvait
@license: GNU General Public License v3.0
"""
import asyncio
import os
import struct
import time
import zipfile
import zlib

import numpy as np
import zarr
from zarr.abc.store import OffsetByteRequest, RangeByteRequest, SuffixByteRequest

_LOCAL_HEADER = struct.Struct("<4s5H3L2H")  # zipfile.structFileHeader, 30 bytes
APPEND_SUFFIX = ".zarzs-append"


def append_record(path):
	"""Name of the file recording an unfinished append to the zip path."""
	return os.fspath(path) + APPEND_SUFFIX


def _append_finished(path, size):
	"""True if the zip has a central directory written after the append that started at size (closed
	without close(), when zipfile wrote it at garbage collection)."""
	try:
		with zipfile.ZipFile(path) as z:
			return z.start_dir >= size
	except (zipfile.BadZipFile, OSError):
		return False


def _check_append_record(path):
	"""Remove the record of an append that finished; raise if one did not finish."""
	record = append_record(path)
	if os.path.exists(record):
		with open(record) as f:
			size = int(f.read().split()[0])
		if not _append_finished(path, size):
			raise zipfile.BadZipFile(f"An append to {path} did not finish ({record} exists); "
									 f"ZARZS.recover(path) restores the zip as it was before")
		os.remove(record)


def recover(path):
	"""Undo an append that did not finish (the process died before close()): truncate the zip to its
	size before the append and remove the record. Returns True if the zip was truncated; an append
	that finished (also without close()) is kept and only its record removed."""
	record = append_record(path)
	if not os.path.exists(record):
		return False
	with open(record) as f:
		size = int(f.read().split()[0])
	if _append_finished(path, size):
		os.remove(record)
		return False
	if os.path.getsize(path) < size:
		raise zipfile.BadZipFile(f"{path} is shorter than before the append recorded in {record}")
	os.truncate(path, size)
	os.remove(record)
	return True


def _fsync(path):
	fd = os.open(path, os.O_RDONLY)
	try:
		os.fsync(fd)
	finally:
		os.close(fd)


class ZipStore(zarr.storage.ZipStore):
	"""zarr.storage.ZipStore that threads write and read in parallel; see the module documentation."""

	supports_deletes = True

	def __init__(self, path, *, mode="r", read_only=None, compression=zipfile.ZIP_STORED, allowZip64=True):
		if mode not in ("r", "w", "a"):
			raise ValueError(f"mode must be 'r', 'w' or 'a', got {mode!r}")
		if mode != "r" and compression != zipfile.ZIP_STORED:
			raise ValueError("ZARZS.ZipStore writes entries uncompressed (compression=zipfile.ZIP_STORED)")
		if mode != "w":
			_check_append_record(path)
		self.replaced_bytes = 0
		super().__init__(path, mode=mode, read_only=read_only, compression=compression, allowZip64=allowZip64)

	# ---- opening, closing, pickling

	def _sync_open(self):
		super()._sync_open()
		self._fd = os.open(self.path, os.O_RDONLY if self._zmode == "r" else os.O_RDWR)
		self._data_offsets = {}
		self._next_offset = 0
		if self._zmode == "a":
			# new entries go after the end of the file: until close() the zip as it was stays intact
			size = os.fstat(self._fd).st_size
			_check_append_record(self.path)
			with open(append_record(self.path), "w") as f:
				f.write(f"{size}\n")
				f.flush()
				os.fsync(f.fileno())
			self.replaced_bytes = size - self._zf.start_dir  # the old central directory
			self._next_offset = size
			self._zf.start_dir = size
		elif self._zmode == "w" and os.path.exists(append_record(self.path)):
			os.remove(append_record(self.path))  # of a zip that is replaced now

	def close(self):
		was_open = self._is_open
		if was_open:
			with self._lock:
				if self._zmode != "r":
					self._zf.start_dir = self._next_offset  # the central directory follows the last entry
				os.close(self._fd)
		super().close()
		if was_open and self._zmode == "a":
			_fsync(self.path)  # the new central directory is on disk before the record goes
			os.remove(append_record(self.path))

	def __getstate__(self):
		if self._zmode != "r":
			raise TypeError("ZARZS.ZipStore open for writing can not be copied to another process "
							"(it would truncate or corrupt the zip); write from threads of one process")
		state = super().__getstate__()
		for name in ("_fd", "_data_offsets", "_next_offset"):
			state.pop(name, None)
		return state

	# ---- writing

	def _set(self, key, value):
		# the bytes of the buffer without a copy (to_bytes would copy the whole chunk holding the GIL)
		data = memoryview(np.ascontiguousarray(value.as_numpy_array())).cast("B")
		zinfo = zipfile.ZipInfo(filename=key, date_time=time.localtime(time.time())[:6])
		zinfo.compress_type = zipfile.ZIP_STORED
		zinfo.external_attr = 0o644 << 16
		zinfo.file_size = zinfo.compress_size = data.nbytes
		zinfo.CRC = zlib.crc32(data)  # releases the GIL for large buffers
		header = zinfo.FileHeader(zip64=data.nbytes > zipfile.ZIP64_LIMIT)
		with self._lock:
			zinfo.header_offset = offset = self._next_offset
			self._next_offset += len(header) + data.nbytes
			# where zipfile writes the central directory, also if it is closed without close() (at exit)
			self._zf.start_dir = self._next_offset
		os.pwrite(self._fd, header, offset)
		written = 0
		while written < data.nbytes:  # pwrite may write less than asked for
			written += os.pwrite(self._fd, data[written:], offset + len(header) + written)
		with self._lock:
			self._remove_entry(key)
			self._zf.filelist.append(zinfo)
			self._zf.NameToInfo[key] = zinfo
			self._data_offsets[key] = offset + len(header)
			self._zf._didModify = True

	def _remove_entry(self, key):
		"""Remove key from the directory of the zip (under the lock); its bytes stay unused in the file."""
		old = self._zf.NameToInfo.pop(key, None)
		if old is not None:
			self._zf.filelist.remove(old)
			self._data_offsets.pop(key, None)
			self.replaced_bytes += old.compress_size
			self._zf._didModify = True
		return old is not None

	async def set(self, key, value):
		self._check_writable()
		if not self._is_open:
			self._sync_open()
		await asyncio.to_thread(self._set, key, value)

	async def set_if_not_exists(self, key, value):
		self._check_writable()
		if not self._is_open:
			self._sync_open()
		with self._lock:
			exists = key in self._zf.NameToInfo
		if not exists:
			await asyncio.to_thread(self._set, key, value)

	async def delete(self, key):
		"""Delete key; if it is no entry, everything under key/ (as zarr's LocalStore deletes a directory)."""
		self._check_writable()
		if not self._is_open:
			self._sync_open()
		with self._lock:
			removed = self._remove_entry(key)
		if not removed:
			await self.delete_dir(key)

	async def delete_dir(self, prefix):
		self._check_writable()
		if not self._is_open:
			self._sync_open()
		if prefix != "" and not prefix.endswith("/"):
			prefix += "/"
		with self._lock:
			for key in [k for k in self._zf.NameToInfo if k.startswith(prefix)]:
				self._remove_entry(key)

	async def clear(self):
		await self.delete_dir("")

	# ---- reading

	def _data_offset(self, zinfo):
		"""Offset of the data of an entry: after its local header (whose name and extra field may differ
		from the central directory)."""
		offset = self._data_offsets.get(zinfo.filename)
		if offset is None:
			header = _LOCAL_HEADER.unpack(os.pread(self._fd, _LOCAL_HEADER.size, zinfo.header_offset))
			if header[0] != zipfile.stringFileHeader:
				raise zipfile.BadZipFile(f"Bad local header of {zinfo.filename} in {self.path}")
			offset = zinfo.header_offset + _LOCAL_HEADER.size + header[-2] + header[-1]
			self._data_offsets[zinfo.filename] = offset
		return offset

	def _get(self, key, prototype, byte_range=None):
		if not self._is_open:
			self._sync_open()
		with self._lock:
			zinfo = self._zf.NameToInfo.get(key)
			if zinfo is not None and zinfo.compress_type != zipfile.ZIP_STORED:
				# a compressed entry of an existing zip: through zipfile
				return super()._get(key, prototype, byte_range)
		if zinfo is None:
			return None
		size = zinfo.file_size
		begin, end = 0, size
		if isinstance(byte_range, RangeByteRequest):
			begin, end = byte_range.start, min(byte_range.end, size)
		elif isinstance(byte_range, OffsetByteRequest):
			begin = byte_range.offset
		elif isinstance(byte_range, SuffixByteRequest):
			begin = max(0, size - byte_range.suffix)
		elif byte_range is not None:
			raise TypeError(f"Unexpected byte_range, got {byte_range}.")
		start = self._data_offset(zinfo)
		return prototype.buffer.from_bytes(os.pread(self._fd, max(0, end - begin), start + begin))

	async def get(self, key, prototype, byte_range=None):
		return await asyncio.to_thread(self._get, key, prototype, byte_range)

	async def get_partial_values(self, prototype, key_ranges):
		return await asyncio.to_thread(lambda: [self._get(key, prototype, byte_range) for key, byte_range in key_ranges])

	async def exists(self, key):
		if not self._is_open:
			self._sync_open()
		with self._lock:
			return key in self._zf.NameToInfo

	async def list(self):
		if not self._is_open:
			self._sync_open()
		with self._lock:
			keys = list(self._zf.NameToInfo)
		for key in keys:
			yield key

	def __repr__(self):
		return f"ZARZS.ZipStore('zip://{self.path}', mode={self._zmode!r})"
