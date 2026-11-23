#!/usr/bin/env python
#-*- coding: utf-8 -*-
"""
Created 2026

Test of ZARZS.ZipStore: parallel writing and reading of zarr zip stores, standard zip output,
replacing and deleting entries, appending to existing zips, zip64.

@author: Vojtech Kulvait
"""
import os
import pickle
import subprocess
import shutil
import tempfile
import unittest
import zipfile
from multiprocessing.dummy import Pool
from unittest import mock

import numpy as np
import zarr

from denpy import ZARZS


def frame(k, shape=(64, 96)):
	return (np.arange(np.prod(shape), dtype=np.uint16).reshape(shape) * (k + 1)).astype(np.uint16)


class TestZipStore(unittest.TestCase):

	def setUp(self):
		self.tmp = tempfile.TemporaryDirectory()
		self.zip = os.path.join(self.tmp.name, "a.zip")

	def tearDown(self):
		self.tmp.cleanup()

	def write(self, n=24, threads=6, path=None, **create):
		store = ZARZS.ZipStore(path or self.zip, mode="w")
		group = zarr.open_group(store, mode="w", attributes={"a": 1})
		params = group.create_group("params", attributes={"b": [1, 2], "c": "text"})
		params.create_array("s", shape=(3,), dtype="string")[:] = np.array(["x", "y", "z"])
		array = group.create_array("img", shape=(n, 64, 96), chunks=(1, 64, 96), dtype=np.uint16, **create)
		with Pool(threads) as pool:
			pool.map(lambda k: array.__setitem__(k, frame(k)), range(n))
		store.close()
		return store

	def check(self, n=24, path=None):
		path = path or self.zip
		with zipfile.ZipFile(path) as z:
			self.assertIsNone(z.testzip())  # CRC-32 of every entry
			names = z.namelist()
		self.assertEqual(len(names), len(set(names)))
		for store in (zarr.storage.ZipStore(path, mode="r"), ZARZS.ZipStore(path, mode="r")):
			group = zarr.open_group(store, mode="r")
			self.assertEqual(group.attrs["a"], 1)
			self.assertEqual(dict(group["params"].attrs), {"b": [1, 2], "c": "text"})
			self.assertEqual(list(group["params/s"][:]), ["x", "y", "z"])
			np.testing.assert_array_equal(group["img"][:], np.stack([frame(k) for k in range(n)]))
			store.close()
		if shutil.which("unzip"):
			out = subprocess.run(["unzip", "-tq", path], capture_output=True, text=True)
			self.assertEqual(out.returncode, 0, out.stdout + out.stderr)

	def test_parallel_write_is_a_standard_zip(self):
		self.write()
		self.check()

	def test_compressed_chunks(self):
		self.write(compressors=zarr.codecs.BloscCodec(cname="zstd", clevel=3, shuffle="bitshuffle"))
		self.check()

	def test_parallel_read(self):
		self.write()
		store = ZARZS.ZipStore(self.zip, mode="r")
		array = zarr.open_array(store, path="img", mode="r")
		with Pool(6) as pool:
			frames = pool.map(lambda k: array[k], range(24))
		np.testing.assert_array_equal(np.stack(frames), np.stack([frame(k) for k in range(24)]))
		# partial reads of a chunk
		from zarr.abc.store import OffsetByteRequest, RangeByteRequest, SuffixByteRequest
		from zarr.core.buffer import default_buffer_prototype
		whole = store._get("img/c/3/0/0", default_buffer_prototype()).to_bytes()
		for request, expected in ((RangeByteRequest(10, 20), whole[10:20]), (OffsetByteRequest(100), whole[100:]),
								  (SuffixByteRequest(7), whole[-7:])):
			self.assertEqual(store._get("img/c/3/0/0", default_buffer_prototype(), request).to_bytes(), expected)
		self.assertIsNone(store._get("missing", default_buffer_prototype()))
		with self.assertRaises(ValueError):  # read only
			zarr.open_group(store, mode="r").attrs.update({"x": 1})
		store.close()

	def test_rewritten_metadata_once(self):
		# group attributes updated after writing: zarr's ZipStore adds a duplicate zarr.json
		store = ZARZS.ZipStore(self.zip, mode="w")
		group = zarr.open_group(store, mode="w")
		group.create_array("x", shape=(4,), dtype=np.int32)[:] = np.arange(4)
		group.attrs.update({"late": True})
		group.attrs.update({"later": 2})
		self.assertGreater(store.replaced_bytes, 0)
		store.close()
		with zipfile.ZipFile(self.zip) as z:
			names = z.namelist()
			self.assertIsNone(z.testzip())
		self.assertEqual(names.count("zarr.json"), 1)
		group = zarr.open_group(zarr.storage.ZipStore(self.zip, mode="r"), mode="r")
		self.assertEqual(dict(group.attrs), {"late": True, "later": 2})

	def test_overwrite_array(self):
		# deleting entries lets zarr overwrite an array (NotImplementedError with zarr's ZipStore)
		store = ZARZS.ZipStore(self.zip, mode="w")
		group = zarr.open_group(store, mode="w")
		group.create_array("x", shape=(8,), chunks=(2,), dtype=np.int32)[:] = np.arange(8)
		group.create_array("x", shape=(3,), chunks=(3,), dtype=np.float32, overwrite=True)[:] = 1.5
		store.close()
		with zipfile.ZipFile(self.zip) as z:
			names = z.namelist()
		self.assertEqual(sorted(names), ["x/c/0", "x/zarr.json", "zarr.json"])
		np.testing.assert_array_equal(zarr.open_array(zarr.storage.ZipStore(self.zip, mode="r"), path="x", mode="r")[:], [1.5] * 3)

	def test_delete_array_by_path(self):
		# delete("x") removes the array x as zarr's LocalStore does (removeHotPixels.py uses it)
		import asyncio
		self.write(n=4)
		store = ZARZS.ZipStore(self.zip, mode="a")
		asyncio.run(store.delete("img"))
		asyncio.run(store.delete("params/s/zarr.json"))  # a single entry
		asyncio.run(store.delete("does/not/exist"))
		store.close()
		with zipfile.ZipFile(self.zip) as z:
			names = z.namelist()
			self.assertIsNone(z.testzip())
		self.assertFalse([n for n in names if n.startswith("img/")])
		self.assertNotIn("params/s/zarr.json", names)
		self.assertIn("params/s/c/0", names)
		self.assertEqual(zarr.open_group(zarr.storage.ZipStore(self.zip, mode="r"), mode="r").attrs["a"], 1)

	def test_without_close(self):
		# a store that is not closed: zipfile writes the central directory when it is garbage
		# collected, after the last entry and not over the entries
		import gc
		for mode in ("w", "a"):
			with self.subTest(mode=mode):
				if mode == "w":
					store = ZARZS.ZipStore(self.zip, mode="w")
					group = zarr.open_group(store, mode="w", attributes={"a": 1})
				else:
					self.write(n=2)
					store = ZARZS.ZipStore(self.zip, mode="a")
					group = zarr.open_group(store, mode="a")
				array = group.create_array("late", shape=(3, 64, 96), chunks=(1, 64, 96), dtype=np.uint16)
				for k in range(3):
					array[k] = frame(k)
				del group, array, store
				gc.collect()
				with zipfile.ZipFile(self.zip) as z:
					self.assertIsNone(z.testzip())
				group = zarr.open_group(zarr.storage.ZipStore(self.zip, mode="r"), mode="r")
				np.testing.assert_array_equal(group["late"][:], np.stack([frame(k) for k in range(3)]))
				if mode == "a":
					np.testing.assert_array_equal(group["img"][:], np.stack([frame(k) for k in range(2)]))
					# the append finished at garbage collection: its record goes when the zip is opened
					ZARZS.ZipStore(self.zip, mode="r")
					self.assertFalse(os.path.exists(ZARZS.append_record(self.zip)))

	def test_append(self):
		self.write()
		store = ZARZS.ZipStore(self.zip, mode="a")
		group = zarr.open_group(store, mode="a")
		self.assertEqual(group.attrs["a"], 1)  # existing entries are read
		extra = group.create_array("extra", shape=(4, 64, 96), chunks=(1, 64, 96), dtype=np.uint16)
		with Pool(4) as pool:
			pool.map(lambda k: extra.__setitem__(k, frame(k + 100)), range(4))
		store.close()
		self.check()
		group = zarr.open_group(zarr.storage.ZipStore(self.zip, mode="r"), mode="r")
		np.testing.assert_array_equal(group["extra"][:], np.stack([frame(k + 100) for k in range(4)]))

	def test_append_keeps_the_old_zip_until_close(self):
		self.write(n=4)
		size = os.path.getsize(self.zip)
		with open(self.zip, "rb") as f:
			original = f.read()
		store = ZARZS.ZipStore(self.zip, mode="a")
		group = zarr.open_group(store, mode="a")
		group.create_array("extra", shape=(3, 64, 96), chunks=(1, 64, 96), dtype=np.uint16)[:] = 5
		self.assertTrue(os.path.exists(ZARZS.append_record(self.zip)))
		with open(self.zip, "rb") as f:
			self.assertEqual(f.read(size), original)  # nothing written over the old zip
		store.close()
		self.assertFalse(os.path.exists(ZARZS.append_record(self.zip)))
		self.check(n=4)

	def test_killed_append_is_recovered(self):
		# a process appending dies before close(): the zip is refused until ZARZS.recover() truncates
		# it, which gives back exactly the zip before the append
		import hashlib
		import sys
		self.write(n=4)
		with open(self.zip, "rb") as f:
			before = hashlib.sha256(f.read()).hexdigest()
		code = (
			"import os, numpy as np, zarr\n"
			"from denpy import ZARZS\n"
			f"store = ZARZS.ZipStore({self.zip!r}, mode='a')\n"
			"group = zarr.open_group(store, mode='a')\n"
			"a = group.create_array('extra', shape=(50, 128, 128), chunks=(1, 128, 128), dtype=np.uint16)\n"
			"a[:] = 7\n"
			"group.attrs.update({'a': 2})\n"
			"os._exit(1)\n")
		root = os.path.dirname(os.path.dirname(os.path.abspath(ZARZS.__file__)))
		env = dict(os.environ, PYTHONPATH=root + os.pathsep + os.environ.get("PYTHONPATH", ""))
		result = subprocess.run([sys.executable, "-c", code], env=env, capture_output=True, text=True)
		self.assertEqual(result.returncode, 1, result.stderr)
		self.assertTrue(os.path.exists(ZARZS.append_record(self.zip)))
		for mode in ("r", "a"):
			with self.assertRaises(zipfile.BadZipFile):
				ZARZS.ZipStore(self.zip, mode=mode)
		self.assertTrue(ZARZS.recover(self.zip))
		with open(self.zip, "rb") as f:
			self.assertEqual(hashlib.sha256(f.read()).hexdigest(), before)
		self.assertFalse(os.path.exists(ZARZS.append_record(self.zip)))
		self.assertFalse(ZARZS.recover(self.zip))
		self.check(n=4)

	def test_compressed_entries_of_an_existing_zip(self):
		# a zip whose entries are deflated (written by another tool) is read through zipfile
		self.write(n=4)
		deflated = os.path.join(self.tmp.name, "deflated.zip")
		with zipfile.ZipFile(self.zip) as src, zipfile.ZipFile(deflated, "w", zipfile.ZIP_DEFLATED) as dst:
			for name in src.namelist():
				dst.writestr(name, src.read(name))
		self.check(n=4, path=deflated)

	def test_zip64(self):
		# central directory beyond the zip64 limit, without writing gigabytes
		with mock.patch.object(zipfile, "ZIP64_LIMIT", 2000):
			self.write(n=8)
		self.check(n=8)

	def test_restrictions(self):
		with self.assertRaises(ValueError):
			ZARZS.ZipStore(self.zip, mode="w", compression=zipfile.ZIP_DEFLATED)
		with self.assertRaises(ValueError):
			ZARZS.ZipStore(self.zip, mode="x")
		store = ZARZS.ZipStore(self.zip, mode="w")
		zarr.open_group(store, mode="w")
		with self.assertRaises(TypeError):  # a writer in another process would truncate the zip
			pickle.dumps(store)
		store.close()
		store = ZARZS.ZipStore(self.zip, mode="r")
		zarr.open_group(store, mode="r")
		copy = pickle.loads(pickle.dumps(store))  # readers can be copied
		self.assertIsNotNone(zarr.open_group(copy, mode="r"))
		store.close()
		copy.close()


if __name__ == "__main__":
	unittest.main()
