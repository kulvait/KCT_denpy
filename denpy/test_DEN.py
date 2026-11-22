#!/usr/bin/env python
#-*- coding: utf-8 -*-
"""
Created 2026

Test of DEN frame access: getFrame and writeFrame for 2D to 5D files.

@author: Vojtech Kulvait
"""
import os
import tempfile
import unittest
import numpy as np

from denpy import DEN


class TestFrameAccess(unittest.TestCase):

	def setUp(self):
		self.tmp = tempfile.TemporaryDirectory()

	def tearDown(self):
		self.tmp.cleanup()

	def store(self, array, name="a.den"):
		path = os.path.join(self.tmp.name, name)
		DEN.storeNdarrayAsDEN(path, array, force=True)
		return path

	def test_flat_index_all_dimensions(self):
		rng = np.random.default_rng(1)
		for shape in ((4, 5, 6), (2, 3, 5, 6), (2, 3, 2, 5, 6)):
			array = rng.random(shape).astype(np.float32)
			path = self.store(array)
			frames = array.reshape((-1,) + shape[-2:])
			for k in range(frames.shape[0]):
				np.testing.assert_array_equal(DEN.getFrame(path, k), frames[k])
				np.testing.assert_array_equal(DEN.getFrame(path, np.int64(k)), frames[k])

	def test_list_index_is_dimspec_order(self):
		array = np.arange(2 * 3 * 4 * 5, dtype=np.float32).reshape(2, 3, 4, 5)  # numpy (z2, z1, y, x)
		path = self.store(array)
		for z2 in range(2):
			for z1 in range(3):
				# list in DEN dimspec order: (z1, z2), the first changes fastest
				np.testing.assert_array_equal(DEN.getFrame(path, [z1, z2]), array[z2, z1])
				np.testing.assert_array_equal(DEN.getFrame(path, (z1, z2)), array[z2, z1])
				# and the same frame by its flat number
				np.testing.assert_array_equal(DEN.getFrame(path, z2 * 3 + z1), array[z2, z1])

	def test_subframe(self):
		array = np.arange(2 * 3 * 4 * 5, dtype=np.float32).reshape(2, 3, 4, 5)
		path = self.store(array)
		np.testing.assert_array_equal(DEN.getFrame(path, 4, row_from=1, row_to=3, col_from=2, col_to=5), array[1, 1, 1:3, 2:5])

	def test_2d(self):
		array = np.arange(20, dtype=np.uint16).reshape(4, 5)
		path = self.store(array)
		np.testing.assert_array_equal(DEN.getFrame(path, 0), array)
		with self.assertRaises(IndexError):
			DEN.getFrame(path, 1)

	def test_invalid_indices(self):
		path4 = self.store(np.zeros((2, 3, 4, 5), dtype=np.float32), "four.den")
		path3 = self.store(np.zeros((3, 4, 5), dtype=np.float32), "three.den")
		with self.assertRaises(IndexError):
			DEN.getFrame(path4, 6)
		with self.assertRaises(IndexError):
			DEN.getFrame(path4, -1)
		with self.assertRaises(IndexError):
			DEN.getFrame(path4, [3, 0])  # z1 < 3 only
		with self.assertRaises(TypeError):
			DEN.getFrame(path4, [1])
		with self.assertRaises(TypeError):
			DEN.getFrame(path4, True)
		with self.assertRaises(TypeError):
			DEN.getFrame(path3, [1])
		with self.assertRaises(TypeError):
			DEN.getFrame(path4, 1.0)

	def test_write_frame(self):
		array = np.zeros((2, 3, 4, 5), dtype=np.float32)
		path = self.store(array)
		frame = np.full((4, 5), 7, dtype=np.float32)
		DEN.writeFrame(path, 4, frame, force=True)              # flat: z2=1, z1=1
		DEN.writeFrame(path, [2, 0], frame * 2, force=True)     # dimspec order: z1=2, z2=0
		data = DEN.getNumpyArray(path)
		np.testing.assert_array_equal(data[1, 1], frame)
		np.testing.assert_array_equal(data[0, 2], frame * 2)
		self.assertEqual(np.count_nonzero(data), 2 * 20)


if __name__ == "__main__":
	unittest.main()
