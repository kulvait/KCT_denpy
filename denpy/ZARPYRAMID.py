#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created 11 2026

Multiscale pyramids of Z,Y,X stacks written with plain zarr, optionally as OME-NGFF (OME-Zarr) images.

How a pyramid is built is described by PyramidSpec (PyramidSpec.describe() says it in words):
	levels   number of levels including the full resolution level "0" (default 5, the same as the
	         default scale factors 2, 4, 8, 16 of ome-zarr-py)
	factor   integer reduction factor between consecutive levels (default 2)
	axes     "yx": every level reduces Y and X, Z is kept (2D reduction, the ome-zarr-py default)
	         "zyx": every level reduces Z, Y and X (3D reduction)
	method   how a block of factor x factor (x factor) samples becomes one sample:
	         "mean" (average), "nearest" (the first sample of the block), "max" or "min"
The size of a level along a reduced axis is the size of the previous level divided by the factor and
rounded down, as in ome-zarr-py; samples at the end of an axis that do not fill a whole block are
dropped. Levels in which a reduced axis would have size 0 are not created.

PyramidWriter takes the frames in Z order and builds all levels while writing, holding at most
`factor` frames per level, so a stack never has to fit in memory. Every level is computed from the
previous level before storing (no lossy codec in between). For a 3D reduction the Y,X blocks of each
frame are reduced first and then `factor` consecutive frames along Z; for the methods offered this is
the same as reducing the 3D blocks at once.

write_pyramid() writes a whole stack, with OME-NGFF metadata if ome=True: version 0.5 for zarr v3
(metadata in attributes["ome"], arrays with dimension_names) and version 0.4 for zarr v2. The
PyramidSpec is stored in the multiscales metadata ("type", "metadata"), so a file records how its
pyramid was made.

@author: Vojtěch Kulvait
@license: GNU General Public License v3.0
"""
import logging
import os
import shutil
from dataclasses import dataclass, asdict

import numpy as np
import zarr

from denpy import ZARZS

from denpy import ZAR

# Create a logger specific to this module
log = logging.getLogger(__name__)
log.setLevel(logging.INFO)
ch = logging.StreamHandler()
ch.setLevel(logging.INFO)
formatter = logging.Formatter('%(asctime)s - %(name)s:%(lineno)d - %(levelname)s : %(message)s', datefmt='%d.%m.%Y %H:%M:%S')
ch.setFormatter(formatter)
log.addHandler(ch)
log.propagate = False

METHODS = ("mean", "nearest", "max", "min")
AXES = ("yx", "zyx")


@dataclass(frozen=True)
class PyramidSpec:
	"""How a pyramid is built: see the module documentation and describe()."""
	levels: int = 5
	factor: int = 2
	axes: str = "yx"
	method: str = "mean"

	def __post_init__(self):
		if not isinstance(self.levels, (int, np.integer)) or self.levels < 1:
			raise ValueError("levels must be an integer >= 1, got %s" % (self.levels,))
		if not isinstance(self.factor, (int, np.integer)) or self.factor < 2:
			raise ValueError("factor must be an integer >= 2, got %s" % (self.factor,))
		if self.axes not in AXES:
			raise ValueError("axes must be one of %s, got %s" % (AXES, self.axes))
		if self.method not in METHODS:
			raise ValueError("method must be one of %s, got %s" % (METHODS, self.method))

	@property
	def reduces_z(self):
		return self.axes == "zyx"

	def level_shapes(self, shape):
		"""Shapes (z, y, x) of the levels for a stack of the given shape (z, y, x).

		Fewer than `levels` shapes are returned when a reduced axis would become 0."""
		z, y, x = (int(s) for s in shape)
		shapes = [(z, y, x)]
		f = self.factor
		while len(shapes) < self.levels:
			z, y, x = (z // f if self.reduces_z else z), y // f, x // f
			if z == 0 or y == 0 or x == 0:
				break
			shapes.append((z, y, x))
		return shapes

	def level_factors(self, count):
		"""Cumulative reduction factors (fz, fy, fx) of the first `count` levels relative to level 0."""
		return [((self.factor**level if self.reduces_z else 1), self.factor**level, self.factor**level) for level in range(count)]

	def describe(self, shape=None):
		"""How the pyramid is built, in words; with a stack shape (z, y, x) also the level shapes."""
		f = self.factor
		if self.reduces_z:
			reduced, block = "Z, Y and X", "%dx%dx%d blocks" % (f, f, f)
		else:
			reduced, block = "Y and X (Z kept)", "%dx%d blocks in each frame" % (f, f)
		how = {"mean": "replaced by their mean", "nearest": "replaced by their first sample",
			   "max": "replaced by their maximum", "min": "replaced by their minimum"}[self.method]
		text = ("%d levels; each level reduces %s of the previous one by %d: %s are %s; "
				"samples at the end of an axis that do not fill a whole block are dropped" % (self.levels, reduced, f, block, how))
		if shape is not None:
			shapes = self.level_shapes(shape)
			text += "; levels: " + ", ".join("%d: %s" % (i, "x".join(str(s) for s in sh)) for i, sh in enumerate(shapes))
			if len(shapes) < self.levels:
				text += " (only %d levels, a smaller one would have an axis of size 0)" % len(shapes)
		return text


def _reduce_yx(frame, factor, method):
	"""Reduce factor x factor blocks of a Y,X frame (trailing rows/columns dropped)."""
	h, w = frame.shape[0] // factor, frame.shape[1] // factor
	if method == "nearest":
		# a copy: the caller may reuse the frame buffer while a 3D reduction still holds this one
		return np.array(frame[:h * factor:factor, :w * factor:factor])
	blocks = frame[:h * factor, :w * factor].reshape(h, factor, w, factor)
	if method == "mean":
		return blocks.mean(axis=(1, 3), dtype=np.float64)
	if method == "max":
		return blocks.max(axis=(1, 3))
	return blocks.min(axis=(1, 3))


def _reduce_z(frames, method):
	"""Reduce `factor` consecutive (already Y,X reduced) frames along Z."""
	if method == "nearest":
		return frames[0]
	stack = np.stack(frames)
	if method == "mean":
		return stack.mean(axis=0, dtype=np.float64)
	if method == "max":
		return stack.max(axis=0)
	return stack.min(axis=0)


def _cast(frame, dtype):
	"""Cast to dtype; floating values going to an integer type are rounded and clipped."""
	dtype = np.dtype(dtype)
	if np.issubdtype(dtype, np.integer) and not np.issubdtype(frame.dtype, np.integer):
		info = np.iinfo(dtype)
		frame = np.clip(np.rint(frame), info.min, info.max)
	return frame.astype(dtype, copy=False)


def default_chunks(level_shape):
	"""Chunks of a level: 1 frame of 256x256 tiles when Y and X exceed 512, else whole frames."""
	_, h, w = level_shape
	return (1, 256, 256) if (h > 512 and w > 512) else (1, h, w)


class PyramidWriter:
	"""Write a Z,Y,X stack frame by frame into the arrays "0", "1", ... of a zarr group.

	group        zarr group (zarr v2 or v3) to create the level arrays in
	shape        (z, y, x) of the full resolution stack
	dtype        dtype of all levels
	spec         PyramidSpec (default PyramidSpec())
	compressors  as returned by ZAR.get_compressor (zarr v3: list of codecs, image codecs such as
	             jpeg2k become the serializer; zarr v2: numcodecs codec or None)
	chunks       None (default_chunks), a tuple for all levels (clipped to the level shape), or a
	             callable level_shape -> chunks
	Frames must be written in Z order with write_frame(); call close() at the end.
	"""

	def __init__(self, group, shape, dtype, spec=None, *, compressors=None, chunks=None):
		self.group = group
		self.spec = spec if spec is not None else PyramidSpec()
		self.dtype = np.dtype(dtype)
		self.shapes = self.spec.level_shapes(shape)
		zarr_format = group.metadata.zarr_format
		self.arrays = []
		for level, level_shape in enumerate(self.shapes):
			if chunks is None:
				level_chunks = default_chunks(level_shape)
			elif callable(chunks):
				level_chunks = tuple(chunks(level_shape))
			else:
				level_chunks = tuple(min(c, s) for c, s in zip(chunks, level_shape))
			extra = {"dimension_names": ("z", "y", "x")} if zarr_format == 3 else {}
			self.arrays.append(group.create_array(str(level), shape=level_shape, dtype=self.dtype,
												  chunks=level_chunks, **ZAR.codec_arguments(compressors), **extra))
		self._written = [0] * len(self.shapes)
		self._buffers = [[] for _ in self.shapes]
		self.min = None
		self.max = None
		self._closed = False

	def write_frame(self, frame):
		"""Write the next frame (Y, X) of level 0 and everything of the other levels it completes."""
		if self._closed:
			raise RuntimeError("PyramidWriter is closed")
		frame = np.asarray(frame)
		if frame.shape != self.shapes[0][1:]:
			raise ValueError("Frame %d has shape %s, expected %s" % (self._written[0], frame.shape, self.shapes[0][1:]))
		if self._written[0] >= self.shapes[0][0]:
			raise ValueError("More frames than the %d of the stack" % self.shapes[0][0])
		if frame.size > 0:
			fmin, fmax = float(np.nanmin(frame)), float(np.nanmax(frame))
			self.min = fmin if self.min is None else min(self.min, fmin)
			self.max = fmax if self.max is None else max(self.max, fmax)
		self._push(0, frame)

	def _push(self, level, frame):
		index = self._written[level]
		if index < self.shapes[level][0]:
			self.arrays[level][index, :, :] = _cast(frame, self.dtype)
		self._written[level] += 1
		if level + 1 >= len(self.shapes):
			return
		reduced = _reduce_yx(frame, self.spec.factor, self.spec.method)
		if not self.spec.reduces_z:
			self._push(level + 1, reduced)
			return
		buffer = self._buffers[level + 1]
		buffer.append(reduced)
		if len(buffer) == self.spec.factor:
			self._push(level + 1, _reduce_z(buffer, self.spec.method))
			buffer.clear()

	def close(self):
		"""Check that the whole stack was written; frames of incomplete Z blocks are dropped by design."""
		if self._written[0] != self.shapes[0][0]:
			raise RuntimeError("Only %d of %d frames written" % (self._written[0], self.shapes[0][0]))
		self._closed = True
		return self


def ome_ngff_attributes(spec, shapes, *, zarr_format=3, name="image", voxel_size=(1.0, 1.0, 1.0), unit=None,
						window=None):
	"""Group attributes of an OME-NGFF image with the levels "0".."n-1" (n = len(shapes)).

	Version 0.5 (under "ome") for zarr_format 3, 0.4 for zarr_format 2. voxel_size is (z, y, x) of
	level 0, the scale of a level is voxel_size times its reduction factors. unit e.g. "micrometer".
	window (min, max) sets the omero display window, if given."""
	axes = []
	for axis in ("z", "y", "x"):
		entry = {"name": axis, "type": "space"}
		if unit is not None:
			entry["unit"] = unit
		axes.append(entry)
	datasets = []
	for level, factors in enumerate(spec.level_factors(len(shapes))):
		scale = [float(v) * f for v, f in zip(voxel_size, factors)]
		datasets.append({"path": str(level), "coordinateTransformations": [{"type": "scale", "scale": scale}]})
	multiscale = {"name": name, "axes": axes, "datasets": datasets,
				  "type": "denpy.ZARPYRAMID." + spec.method,
				  "metadata": dict(asdict(spec), description=spec.describe())}
	omero = None
	if window is not None:
		vmin, vmax = float(window[0]), float(window[1])
		omero = {"channels": [{"color": "FFFFFF", "label": name, "active": True,
							   "window": {"start": vmin, "end": vmax, "min": vmin, "max": vmax}}]}
	if zarr_format == 3:
		ome = {"version": "0.5", "multiscales": [multiscale]}
		if omero is not None:
			ome["omero"] = omero
		return {"ome": ome}
	attributes = {"multiscales": [dict(multiscale, version="0.4")]}
	if omero is not None:
		attributes["omero"] = dict(omero, version="0.4")
	return attributes


def write_pyramid(output, frames, shape, dtype, spec=None, *, zarr_format=3, compressors=None, chunks=None,
				  overwrite=False, ome=True, name="image", voxel_size=(1.0, 1.0, 1.0), unit=None, verbose=False):
	"""Write a Z,Y,X stack as a pyramid, by default an OME-NGFF image; returns the PyramidWriter.

	output   path of a directory store, a path ending with .zip (zip store, written directly with
	         ZARZS.ZipStore, which keeps every entry once although the group metadata is written again
	         at the end), or a zarr store object
	frames   the frames (Y, X) in Z order: an iterable, or a callable k -> frame
	shape    (z, y, x), dtype, spec, compressors, chunks: see PyramidWriter
	ome      write OME-NGFF metadata (ome_ngff_attributes) with name, voxel_size, unit and the
	         data range as display window
	"""
	spec = spec if spec is not None else PyramidSpec()
	zip_store = None
	if isinstance(output, (str, os.PathLike)):
		path = os.fspath(output)
		if os.path.exists(path):
			if not overwrite:
				raise IOError("%s exists, use overwrite=True" % path)
			if os.path.isdir(path):
				shutil.rmtree(path)
			else:
				os.remove(path)
		if path.endswith(".zip"):
			# ZARZS.ZipStore replaces the group metadata written again at the end (zarr's ZipStore
			# would add a second zarr.json) and lets threads write chunks in parallel
			store = zip_store = ZARZS.ZipStore(path, mode="w")
		else:
			store = path
	else:
		store = output
	if verbose:
		log.info("Pyramid of a stack %s, %s, zarr v%d: %s" % ("x".join(str(s) for s in shape), np.dtype(dtype), zarr_format, spec.describe(shape)))
	group = zarr.open_group(store=store, mode="w", zarr_format=zarr_format)
	writer = PyramidWriter(group, shape, dtype, spec, compressors=compressors, chunks=chunks)
	n = int(shape[0])
	source = (frames(k) for k in range(n)) if callable(frames) else iter(frames)
	for k, frame in enumerate(source):
		writer.write_frame(frame)
		if verbose and (k % 100 == 0 or k == n - 1):
			log.info("Written frame %d/%d" % (k + 1, n))
	writer.close()
	if ome:
		window = (writer.min, writer.max) if writer.min is not None else None
		group.attrs.update(ome_ngff_attributes(spec, writer.shapes, zarr_format=zarr_format, name=name,
											   voxel_size=voxel_size, unit=unit, window=window))
	if zip_store is not None:
		zip_store.close()
	return writer
