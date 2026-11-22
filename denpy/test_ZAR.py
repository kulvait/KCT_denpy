#!/usr/bin/env python
#-*- coding: utf-8 -*-
"""
Created 2026

Test of ZAR.get_compressor: every codec name round-trips through a real zarr v2 and v3 array.

@author: Vojtech Kulvait
"""
import os
import tempfile
import unittest
import numpy as np
import zarr

from denpy import ZAR

GENERAL_V2 = ("none", "zstd", "lz4", "gzip", "blosc", "blosc-blosclz", "blosc-zstd", "blosc-lz4", "blosc-zlib")
GENERAL_V3 = ("none", "zstd", "lz4", "gzip", "blosc", "blosc-blosclz", "blosc-zstd", "blosc-lz4", "blosc-lz4hc",
			  "blosc-snappy", "blosc-zlib")
# image codecs: (name, clevel, dtype, lossless)
IMAGE = (("jpeg2k", 0, np.uint16, True), ("htj2k", 0, np.uint16, True), ("jpegxl", 0, np.uint16, True),
		 ("zfp", 0, np.float32, True), ("jpeg2k", 5, np.uint16, False), ("htj2k", 0.05, np.uint16, False),
		 ("jpegxr", 1, np.uint8, False), ("avif", 100, np.uint8, False), ("sz3", 0.01, np.float32, False))


def smooth(shape, dtype):
	"""A smooth image, so that lossy codecs stay close."""
	y, x = np.mgrid[0:shape[0], 0:shape[1]]
	image = 0.5 + 0.25 * np.sin(x / 7.0) + 0.25 * np.cos(y / 5.0)
	if np.issubdtype(dtype, np.integer):
		return (image * (np.iinfo(dtype).max // 2)).astype(dtype)
	return image.astype(dtype)


class TestGetCompressor(unittest.TestCase):

	def setUp(self):
		self.tmp = tempfile.TemporaryDirectory()

	def tearDown(self):
		self.tmp.cleanup()

	def roundtrip(self, name, clevel, data, zarrv2, chunks):
		try:
			codec = ZAR.get_compressor(name, clevel=clevel, zarrv2=zarrv2, dtype=data.dtype)
		except RuntimeError as e:  # codec not built into this imagecodecs
			self.skipTest(str(e))
		path = os.path.join(self.tmp.name, "%s_%s_v%d.zarr" % (name, clevel, 2 if zarrv2 else 3))
		array = zarr.create_array(store=path, shape=data.shape, chunks=chunks, dtype=data.dtype,
								  zarr_format=2 if zarrv2 else 3, **ZAR.codec_arguments(codec))
		array[...] = data
		return zarr.open_array(path, mode="r")[...]

	def test_general_codecs_lossless(self):
		data = np.arange(4 * 32 * 48, dtype=np.float32).reshape(4, 32, 48)
		for zarrv2, names in ((True, GENERAL_V2), (False, GENERAL_V3)):
			for name in names:
				with self.subTest(name=name, zarrv2=zarrv2):
					np.testing.assert_array_equal(self.roundtrip(name, 5, data, zarrv2, (1, 32, 48)), data)

	def test_image_codecs(self):
		for zarrv2 in (True, False):
			for name, clevel, dtype, lossless in IMAGE:
				with self.subTest(name=name, clevel=clevel, zarrv2=zarrv2):
					data = smooth((64, 96), dtype)
					back = self.roundtrip(name, clevel, data, zarrv2, (64, 96))
					self.assertEqual(back.shape, data.shape)
					self.assertEqual(back.dtype, data.dtype)
					if lossless:
						np.testing.assert_array_equal(back, data)
					else:
						# the quality meaning of clevel differs per codec; only check a sane round trip
						self.assertTrue(np.all(np.isfinite(back.astype(np.float64))))

	def test_blosc_typesize(self):
		for dtype in (np.uint8, np.uint16, np.float32, np.float64):
			itemsize = np.dtype(dtype).itemsize
			v2 = ZAR.get_compressor("blosc-zstd", clevel=5, zarrv2=True, dtype=dtype)
			self.assertEqual(v2._typesize, itemsize)  # numcodecs keeps it private
			v3 = ZAR.get_compressor("blosc-zstd", clevel=5, zarrv2=False, dtype=dtype)[0]
			self.assertEqual(v3.typesize, itemsize)

	def test_read_in_new_process(self):
		"""An array written with an image codec can be read by another process after register_codecs()."""
		import subprocess, sys
		for zarrv2 in (True, False):
			with self.subTest(zarrv2=zarrv2):
				data = smooth((64, 96), np.uint16)
				self.roundtrip("htj2k", 0, data, zarrv2, (64, 96))
				path = os.path.join(self.tmp.name, "htj2k_0_v%d.zarr" % (2 if zarrv2 else 3))
				code = ("import sys, numpy as np, zarr; from denpy import ZAR; ZAR.register_codecs(); "
						"a = zarr.open_array(sys.argv[1], mode='r')[...]; print(int(a.sum()))")
				out = subprocess.run([sys.executable, "-c", code, path], capture_output=True, text=True)
				self.assertEqual(out.returncode, 0, out.stderr)
				self.assertEqual(int(out.stdout.strip()), int(data.sum()))

	def test_codec_arguments(self):
		self.assertEqual(ZAR.codec_arguments(None), {"compressors": None})
		v2 = ZAR.get_compressor("zstd", zarrv2=True, dtype=np.uint16)
		self.assertIs(ZAR.codec_arguments(v2)["compressors"], v2)
		blosc = ZAR.get_compressor("blosc-zstd", zarrv2=False, dtype=np.uint16)
		self.assertEqual(ZAR.codec_arguments(blosc), {"compressors": blosc})
		jp2 = ZAR.get_compressor("jpeg2k", clevel=0, zarrv2=False, dtype=np.uint16)
		self.assertEqual(ZAR.codec_arguments(jp2), {"compressors": None, "serializer": jp2[0]})
		self.assertEqual(ZAR.codec_arguments([]), {"compressors": None})

	def test_unavailable_blosc_compressor(self):
		from numcodecs import blosc
		if "snappy" in blosc.list_compressors():
			self.skipTest("this Blosc build has snappy")
		with self.assertRaises(RuntimeError):
			ZAR.get_compressor("blosc-snappy", zarrv2=False)

	def test_unknown_name(self):
		for zarrv2 in (True, False):
			with self.assertRaises(ValueError):
				ZAR.get_compressor("nonexistent", zarrv2=zarrv2)


if __name__ == "__main__":
	unittest.main()
