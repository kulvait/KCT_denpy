#!/usr/bin/env python
#-*- coding: utf-8 -*-
"""
Created 2026

Test of ZARPYRAMID: pyramid levels against a reference reduction of the whole array, metadata.

@author: Vojtech Kulvait
"""
import os
import tempfile
import unittest
import numpy as np
import zarr

from denpy import ZARPYRAMID
from denpy.ZARPYRAMID import PyramidSpec, write_pyramid


def reference_level(volume, spec, level):
	"""Level `level` computed at once from the whole volume (independent of the streaming writer)."""
	f = spec.factor**level
	fz = f if spec.reduces_z else 1
	z, y, x = volume.shape
	nz, ny, nx = z // fz, y // f, x // f
	v = volume[:nz * fz, :ny * f, :nx * f].reshape(nz, fz, ny, f, nx, f)
	if spec.method == "nearest":
		return volume[:nz * fz:fz, :ny * f:f, :nx * f:f]
	if spec.method == "mean":
		return v.mean(axis=(1, 3, 5), dtype=np.float64)
	if spec.method == "max":
		return v.max(axis=(1, 3, 5))
	return v.min(axis=(1, 3, 5))


class TestPyramid(unittest.TestCase):

	def setUp(self):
		self.tmp = tempfile.TemporaryDirectory()

	def tearDown(self):
		self.tmp.cleanup()

	def path(self, name):
		return os.path.join(self.tmp.name, name)

	def check(self, volume, spec, zarr_format=3, name="p.zarr"):
		write_pyramid(self.path(name), volume, volume.shape, volume.dtype, spec, zarr_format=zarr_format)
		group = zarr.open_group(self.path(name), mode="r")
		shapes = spec.level_shapes(volume.shape)
		for level, shape in enumerate(shapes):
			data = group[str(level)][:]
			self.assertEqual(data.shape, shape)
			expected = reference_level(volume.astype(np.float64) if spec.method == "mean" else volume, spec, level)
			if np.issubdtype(volume.dtype, np.integer) and spec.method == "mean":
				# the streaming mean is rounded at every level, so allow half a unit per level
				np.testing.assert_allclose(data, expected, atol=0.5 * level + 1e-9)
			else:
				np.testing.assert_allclose(data, expected, rtol=1e-6)
		return group

	def test_methods_2d_and_3d(self):
		rng = np.random.default_rng(3)
		volume = rng.random((9, 37, 50)).astype(np.float32)  # odd sizes: trailing samples dropped
		for axes in ("yx", "zyx"):
			for method in ("mean", "nearest", "max", "min"):
				with self.subTest(axes=axes, method=method):
					self.check(volume, PyramidSpec(levels=4, axes=axes, method=method), name="%s_%s.zarr" % (axes, method))

	def test_integer_and_factor_3(self):
		rng = np.random.default_rng(4)
		volume = (rng.random((10, 60, 61)) * 60000).astype(np.uint16)
		self.check(volume, PyramidSpec(levels=3, factor=3, axes="zyx", method="mean"))

	def test_zarr_v2(self):
		volume = np.arange(4 * 16 * 16, dtype=np.float32).reshape(4, 16, 16)
		group = self.check(volume, PyramidSpec(levels=3), zarr_format=2)
		self.assertEqual(group.metadata.zarr_format, 2)
		self.assertEqual(group.attrs["multiscales"][0]["version"], "0.4")

	def test_level_shapes_capped(self):
		spec = PyramidSpec(levels=10, axes="zyx")
		self.assertEqual(spec.level_shapes((5, 64, 64)), [(5, 64, 64), (2, 32, 32), (1, 16, 16)])
		self.assertIn("only 3 levels", spec.describe((5, 64, 64)))
		self.assertEqual(PyramidSpec().level_shapes((5, 64, 64))[-1], (5, 4, 4))  # default: 5 levels, Z kept

	def test_metadata_and_callable_frames(self):
		volume = np.arange(6 * 20 * 30, dtype=np.float32).reshape(6, 20, 30)
		spec = PyramidSpec(levels=3, axes="zyx", method="max")
		write_pyramid(self.path("m.zarr"), lambda k: volume[k], volume.shape, np.float32, spec,
					  name="vol", voxel_size=(2.0, 0.5, 0.5), unit="micrometer")
		group = zarr.open_group(self.path("m.zarr"), mode="r")
		ome = group.attrs["ome"]
		self.assertEqual(ome["version"], "0.5")
		ms = ome["multiscales"][0]
		self.assertEqual([d["path"] for d in ms["datasets"]], ["0", "1", "2"])
		self.assertEqual(ms["datasets"][2]["coordinateTransformations"][0]["scale"], [8.0, 2.0, 2.0])
		self.assertEqual(ms["axes"][0], {"name": "z", "type": "space", "unit": "micrometer"})
		self.assertEqual(ms["metadata"]["method"], "max")
		self.assertEqual(ome["omero"]["channels"][0]["window"]["max"], float(volume.max()))
		self.assertEqual(tuple(group["1"].metadata.dimension_names), ("z", "y", "x"))

	def test_zip_store(self):
		volume = np.ones((2, 8, 8), dtype=np.uint8)
		write_pyramid(self.path("z.zip"), volume, volume.shape, volume.dtype, PyramidSpec(levels=2))
		group = zarr.open_group(zarr.storage.ZipStore(self.path("z.zip"), mode="r"), mode="r")
		self.assertEqual(group["1"].shape, (2, 4, 4))
		self.assertIn("ome", group.attrs)
		import zipfile
		with zipfile.ZipFile(self.path("z.zip")) as archive:
			names = archive.namelist()
			self.assertIsNone(archive.testzip())  # CRC-32 of every entry
		self.assertEqual(len(names), len(set(names)))  # every entry once, also the rewritten zarr.json
		self.assertFalse(os.path.exists(self.path("z.zip.tmpdir")))

	def test_reused_frame_buffer(self):
		volume = np.arange(4 * 8 * 8, dtype=np.float32).reshape(4, 8, 8)
		buffer = np.empty((8, 8), dtype=np.float32)

		def frames():
			for k in range(4):
				buffer[:] = volume[k]  # the same array for every frame
				yield buffer
		spec = PyramidSpec(levels=2, axes="zyx", method="nearest")
		write_pyramid(self.path("b.zarr"), frames(), volume.shape, volume.dtype, spec)
		np.testing.assert_array_equal(zarr.open_group(self.path("b.zarr"), mode="r")["1"][:], reference_level(volume, spec, 1))

	def test_image_codecs(self):
		from denpy import ZAR
		rng = np.random.default_rng(5)
		volume = (rng.random((3, 64, 96)) * 4000).astype(np.uint16)
		for name, zarr_format in (("jpeg2k", 3), ("htj2k", 3), ("jpeg2k", 2), ("htj2k", 2)):
			with self.subTest(codec=name, zarr_format=zarr_format):
				codec = ZAR.get_compressor(name, clevel=0, zarrv2=zarr_format == 2, dtype=volume.dtype)  # lossless
				path = self.path("%s_v%d.zarr" % (name, zarr_format))
				write_pyramid(path, volume, volume.shape, volume.dtype, PyramidSpec(levels=3),
							  zarr_format=zarr_format, compressors=codec)
				group = zarr.open_group(path, mode="r")
				np.testing.assert_array_equal(group["0"][:], volume)
				self.assertEqual(group["2"].shape, (3, 16, 24))

	def test_invalid_spec(self):
		for kwargs in ({"levels": 0}, {"factor": 1}, {"axes": "xy"}, {"method": "median"}):
			with self.assertRaises(ValueError):
				PyramidSpec(**kwargs)


if __name__ == "__main__":
	unittest.main()
