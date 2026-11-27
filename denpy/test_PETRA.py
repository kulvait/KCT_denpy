#!/usr/bin/env python
#-*- coding: utf-8 -*-
"""
Created 2026

Test of PETRA: reading P05 scan HDF5 files in the two layouts of the beamline files and their
border cases (missing or empty beam current, empty or scalar camera values, incomplete logs).

  2021 layout: group attributes are 1-element object arrays, strings are (1,) datasets
  2025 layout: group attributes are plain strings, strings are scalar datasets, some values are
               empty (0,) datasets

@author: Vojtech Kulvait
"""
import os
import pathlib
import tempfile
import unittest

import h5py
import numpy as np
import pandas as pd

from denpy import PETRA

KEYS = [2, 2, 1, 1, 1, 0, 0, 0, 0, 1]
T0 = 1_758_209_568_000


def make_scan(path, style="2025", current=True, mutate=None):
	"""Write a small scan file; mutate(h5) changes it before it is closed. Returns the path."""
	old = style == "2021"
	with h5py.File(path, "w") as h5:
		def attr(obj, value):
			if old:
				obj.attrs.create("NX_class", np.array([value], dtype=object), dtype=h5py.string_dtype())
			else:
				obj.attrs["NX_class"] = value

		def group(name, nx_class="NXcollection"):
			g = h5.require_group(name)
			attr(g, nx_class)
			return g

		def text(name, value):
			if old:
				h5.create_dataset(name, data=np.array([value], dtype=object), dtype=h5py.string_dtype())
			else:
				h5.create_dataset(name, data=value.encode(), dtype=h5py.string_dtype("ascii"))

		group("entry", "NXentry")
		group("entry/beamline")
		text("entry/beamline/experiment", "microtomography")
		text("entry/beamline/name", "P05 IBL EH2")
		camera = group("entry/hardware/camera")
		for name, value in (("pixelsize", 6.5), ("magnification", 10.0), ("roi_width", 100.0), ("roi_height", 50.0)):
			camera.create_dataset(name, data=np.array([value]))
		setup = group("entry/scan/setup")
		setup.create_dataset("pos_o_ccd_dist", data=np.array([200.0]))
		setup.create_dataset("pos_p05_energy", data=np.array([20000.0], dtype=np.float32))
		data = group("entry/scan/data")
		times = T0 + np.arange(len(KEYS), dtype=np.uint64) * 100
		files = ["scan_%06d_%s.tif" % (i, {0: "img", 1: "ref", 2: "dar"}[k]) for i, k in enumerate(KEYS)]
		data.create_dataset("image_key/value", data=np.array(KEYS, dtype=np.int16))
		data.create_dataset("image_key/time", data=times)
		data.create_dataset("image_file/value", data=np.array([f.encode() for f in files], dtype=object),
							dtype=h5py.string_dtype("utf-8" if old else "ascii"))
		data.create_dataset("image_file/time", data=times)
		moving = np.array(KEYS) != 2
		data.create_dataset("s_rot/value", data=np.linspace(0, 180, moving.sum()).astype(np.float32))
		data.create_dataset("s_rot/time", data=times[moving])
		data.create_dataset("s_stage_x/value", data=np.linspace(1.0, 2.0, moving.sum()).astype(np.float32))
		data.create_dataset("s_stage_x/time", data=times[moving])
		for name in ("image_key", "image_file", "s_rot", "s_stage_x"):
			attr(data[name], "NXlog")
		if current:
			log = group("entry/hardware/beam_current/current", "NXlog")
			# readings every 250 ms from 100 ms before the scan, 100 mA falling by 1 mA/s
			current_times = T0 - 100 + np.arange(6, dtype=np.uint64) * 250
			log.create_dataset("time", data=current_times)
			log.create_dataset("value", data=100.0 - (current_times.astype(np.float64) - T0) / 1000.0)
		if mutate is not None:
			mutate(h5)
	return path


def replace(h5, path, data, dtype=None):
	del h5[path]
	h5.create_dataset(path, data=data, dtype=dtype)


class TestTimestamps(unittest.TestCase):
	"""fix_timestamps: frames keep their order, only frames with a wrong time are re-timed; known issues
	are recognised (and fixed with fixKnownTimestampIssues), others get the default correction."""

	# the block pattern of P05/P07 files: the flats after the darks and after the projections carry the
	# start time of the block before them
	DAR, REF, IMG, REF2 = np.arange(0, 500, 100), np.arange(0, 300, 100), np.arange(2000, 3000, 100), np.arange(2000, 2300, 100)
	KEYS = np.array([2] * 5 + [1] * 3 + [0] * 10 + [1] * 3)

	def blocks(self):
		return np.concatenate([self.DAR, self.REF, self.IMG, self.REF2])

	def test_flat_field_blocks(self):
		new, issues = PETRA.fix_timestamps(self.blocks(), self.KEYS)
		np.testing.assert_array_equal(new[:5], self.DAR)                        # unchanged
		np.testing.assert_array_equal(new[5:8], [1100, 1200, 1300])            # centred: 700 ms before and after
		np.testing.assert_array_equal(new[8:18], self.IMG)                     # unchanged: the times fit again
		np.testing.assert_array_equal(new[18:], [3600, 3700, 3800])            # the pause before the projections
		self.assertEqual([(k, r) for k, r, _ in issues], [("flat_blocks", (5, 7)), ("flat_blocks", (18, 20))])
		messages = PETRA.timestamp_issue_messages(issues)
		self.assertEqual(len(messages), 1)
		self.assertTrue(messages[0].startswith(PETRA.TIMESTAMP_FLAT_BLOCKS))
		self.assertIn("rows 5-7 +1100 ms, rows 18-20 +1600 ms; fixed", messages[0])
		# without the fix: reported, default correction (1 ms behind the frame before, steps kept)
		new, issues = PETRA.fix_timestamps(self.blocks(), self.KEYS, fixKnownTimestampIssues=False)
		np.testing.assert_array_equal(new[5:8], [401, 501, 601])
		np.testing.assert_array_equal(new[18:], [2901, 3001, 3101])
		self.assertIn("default correction (fixKnownTimestampIssues=False)", PETRA.timestamp_issue_messages(issues, False)[0])
		# without image keys the structure is unknown: default correction
		new, issues = PETRA.fix_timestamps(self.blocks())
		np.testing.assert_array_equal(new[5:8], [401, 501, 601])
		self.assertEqual({k for k, _, _ in issues}, {"unknown"})
		self.assertTrue(PETRA.timestamp_issue_messages(issues)[0].startswith(PETRA.TIMESTAMP_UNKNOWN))

	def test_frame_with_the_time_of_an_earlier_frame(self):
		# a projection with the time of the frame 20 rows earlier, its neighbours on time
		old = np.arange(0, 4000, 100)
		old[30] = old[10]
		new, issues = PETRA.fix_timestamps(old, np.zeros(40))
		np.testing.assert_array_equal(new, np.arange(0, 4000, 100))          # back on the line
		self.assertEqual(issues, [("repeated 20", (30, 30), "row 30 +2000 ms (time of row 10)")])
		self.assertEqual(PETRA.timestamp_issue_messages(issues),
						 [PETRA.TIMESTAMP_REPEATED % 20 + ": row 30 +2000 ms (time of row 10); fixed"])
		new, _ = PETRA.fix_timestamps(old, fixKnownTimestampIssues=False)
		self.assertEqual(new[30], 2901)                                         # default: 1 ms behind
		# a repeated neighbour is the same issue with offset 1
		new, issues = PETRA.fix_timestamps([0, 100, 100, 200, 300])
		np.testing.assert_array_equal(new, [0, 100, 150, 200, 300])
		self.assertEqual(issues[0][0], "repeated 1")

	def test_unknown(self):
		for old, expected in (([0, 100, 100, 100, 101, 500], [0, 100, 101, 102, 103, 500]),   # several frames
							  ([0, 100, 200, 50], [0, 100, 200, 201]),                        # the last frame
							  ([0, 100, 200, 150, 160, 400], [0, 100, 200, 201, 211, 400])):     # a drop of two
			with self.subTest(old=old):
				new, issues = PETRA.fix_timestamps(old)
				np.testing.assert_array_equal(new, expected)
				self.assertEqual({k for k, _, _ in issues}, {"unknown"})
		self.assertEqual(len(PETRA.fix_timestamps([])[0]), 0)
		np.testing.assert_array_equal(PETRA.normalize_timestamps_preserve_increase([7]), [7])

	def test_scan_dataset(self):
		# the flats after the darks logged with the darks' start, the final flat with the projections' start
		times = T0 + np.array([0, 100, 0, 100, 200, 1000, 1100, 1200, 1300, 1000], dtype=np.uint64)

		def mutate(h5):
			for name in ("image_key", "image_file"):
				replace(h5, "entry/scan/data/%s/time" % name, times)
			moving = np.array(KEYS) != 2
			for name in ("s_rot", "s_stage_x"):
				replace(h5, "entry/scan/data/%s/time" % name, times[moving])
		with tempfile.TemporaryDirectory() as tmp:
			path = make_scan(os.path.join(tmp, "scan.h5"), mutate=mutate)
			with self.assertLogs(PETRA.log, "WARNING") as logged:
				df = PETRA.scanDataset(path, includeCurrent=True)
			with h5py.File(path, "r") as h5:
				names = [n.decode() for n in h5["entry/scan/data/image_file/value"][()]]
			self.assertEqual(list(df["image_file"]), names)                    # the order of the HDF5 file
			# flats centred between darks and projections (gaps 350 ms), final flat 350 ms after the projections
			np.testing.assert_array_equal(df.index.to_numpy().astype(np.int64) - T0, [0, 100, 450, 550, 650, 1000, 1100, 1200, 1300, 1650])
			self.assertTrue(df["s_rot"].iloc[2:].notna().all())                # logs at the original times are inserted
			self.assertTrue(any(PETRA.TIMESTAMP_FLAT_BLOCKS in m for m in logged.output))
			df = PETRA.scanDataset(path, fixKnownTimestampIssues=False)
			np.testing.assert_array_equal(df.index.to_numpy().astype(np.int64) - T0, [0, 100, 101, 201, 301, 1000, 1100, 1200, 1300, 1301])
			# getExperimentInfo: end time and duration from the corrected times
			info = PETRA.getExperimentInfo(path)
			self.assertAlmostEqual(info["scan_duration_sec"], 1.65)
			self.assertAlmostEqual(info["sample_duration_sec"], 0.3)
			self.assertAlmostEqual(PETRA.getExperimentInfo(path, fixKnownTimestampIssues=False)["scan_duration_sec"], 1.301)


class TestPETRA(unittest.TestCase):

	def setUp(self):
		self.tmp = tempfile.TemporaryDirectory()

	def tearDown(self):
		self.tmp.cleanup()

	def scan(self, name="scan", **kwargs):
		return make_scan(os.path.join(self.tmp.name, name + ".h5"), **kwargs)

	def test_both_layouts(self):
		for style in ("2021", "2025"):
			with self.subTest(style=style):
				path = self.scan(style=style)
				info = PETRA.getExperimentInfo(path)
				self.assertEqual((info["experiment"], info["beamline"]), ("microtomography", "P05 IBL EH2"))
				self.assertAlmostEqual(info["pix_size"], 0.65)
				self.assertAlmostEqual(info["field_of_view_x"], 65.0)
				self.assertEqual((info["energy_keV"], info["propagation_distance_mm"]), (20.0, 200.0))
				self.assertEqual(info["sample_frame_count"], 4)
				df = PETRA.scanDataset(path)
				self.assertEqual(list(df["image_key"]), KEYS)
				self.assertEqual(df["image_file"].iloc[0], "scan_000000_dar.tif")
				self.assertTrue(df["s_rot"].iloc[:2].isna().all() and df["s_rot"].iloc[2:].notna().all())
				images = PETRA.imageDataset(path, includeCurrent=False, includePixelShift=True)
				self.assertEqual(len(images), 4)
				self.assertIn("pixel_shift", images.columns)

	def test_current_interpolation(self):
		df = PETRA.scanDataset(self.scan(), includeCurrent=True)
		# linear between readings, the readings fall by 1 mA/s
		expected = 100.0 - (df.index.to_numpy().astype(np.float64) - T0) / 1000.0
		np.testing.assert_allclose(df["current"].to_numpy(), expected)
		# a time offset shifts the readings
		shifted = PETRA.scanDataset(self.scan(), includeCurrent=True, timeOffsetSec=0.5)
		# the readings start 400 ms after the scan now: the first 4 frames get the first reading
		np.testing.assert_allclose(shifted["current"].to_numpy()[4:], expected[4:] + 0.5)
		np.testing.assert_allclose(shifted["current"].to_numpy()[:4], 100.1)
		# a pathlib.Path works as a file name
		PETRA.scanDataset(pathlib.Path(self.scan()), includeCurrent=True)

	def test_current_outside_readings_and_nan(self):
		def few(h5):  # one reading in the middle of the scan, one NaN reading
			replace(h5, "entry/hardware/beam_current/current/time", np.array([T0 + 400, T0 + 500], dtype=np.uint64))
			replace(h5, "entry/hardware/beam_current/current/value", np.array([90.0, np.nan]))
		df = PETRA.scanDataset(self.scan(mutate=few), includeCurrent=True)
		np.testing.assert_allclose(df["current"].to_numpy(), 90.0)
		self.assertEqual(len(PETRA.beamCurrentDataset(self.scan(mutate=few))), 1)

	def test_missing_current_is_a_clean_error(self):
		def empty(h5):
			replace(h5, "entry/hardware/beam_current/current/time", np.array([], dtype=np.uint64))
			replace(h5, "entry/hardware/beam_current/current/value", np.array([], dtype=np.float64))

		def all_nan(h5):
			replace(h5, "entry/hardware/beam_current/current/value", np.full(6, np.nan))

		def mismatch(h5):
			replace(h5, "entry/hardware/beam_current/current/value", np.zeros(5))
		cases = {"no current": dict(current=False), "empty": dict(mutate=empty), "all NaN": dict(mutate=all_nan),
				 "length mismatch": dict(mutate=mismatch)}
		for case, kwargs in cases.items():
			with self.subTest(case=case):
				path = self.scan(**kwargs)
				with self.assertRaises(PETRA.BeamCurrentError):
					PETRA.scanDataset(path, includeCurrent=True)
				with self.assertRaises(PETRA.BeamCurrentError):
					PETRA.imageDataset(path)
				# without the current the file is fine
				self.assertEqual(len(PETRA.scanDataset(path)), len(KEYS))
		self.assertTrue(issubclass(PETRA.BeamCurrentError, ValueError))
		self.assertEqual(PETRA.findInsertionEvents(self.scan(mutate=empty)), [])
		with self.assertRaises(PETRA.BeamCurrentError):
			PETRA.findInsertionEvents(self.scan(current=False))

	def test_camera_values(self):
		def empty_magnification(h5):
			replace(h5, "entry/hardware/camera/magnification", np.array([], dtype=np.float64))

		def zero_magnification(h5):
			replace(h5, "entry/hardware/camera/magnification", np.array([0.0]))

		def scalars(h5):
			replace(h5, "entry/hardware/camera/pixelsize", 6.5)
			replace(h5, "entry/hardware/camera/magnification", 5.0)

		def calibration(h5):
			replace(h5, "entry/hardware/camera/magnification", np.array([], dtype=np.float64))
			h5.create_dataset("entry/hardware/camera/calibration/magnification", data=np.array([20.0]))

		def old_names(h5):
			del h5["entry/hardware/camera"]
			camera = h5.create_group("entry/hardware/camera1")
			camera.create_dataset("px_size", data=np.array([6.5]))
			camera.create_dataset("magnification", data=np.array([10.0]))
			camera.create_dataset("senzor_xsize", data=np.array([200.0]))
		for style in ("2021", "2025"):
			for mutate in (empty_magnification, zero_magnification):
				with self.subTest(style=style, case=mutate.__name__):
					path = self.scan(style=style, mutate=mutate)
					info = PETRA.getExperimentInfo(path)
					self.assertNotIn("pix_size", info)
					self.assertNotIn("fresnel_number", info)
					with self.assertRaises(ValueError):
						PETRA.imageDataset(path, includeCurrent=False, includePixelShift=True)
					# an explicit magnification is used instead
					self.assertAlmostEqual(PETRA.getExperimentInfo(path, overrideMagnification=13.0)["pix_size"], 0.5)
					PETRA.imageDataset(path, includeCurrent=False, includePixelShift=True, overrideMagnification=13.0)
			for mutate, pix_size in ((scalars, 1.3), (calibration, 0.325), (old_names, 0.65)):
				with self.subTest(style=style, case=mutate.__name__):
					path = self.scan(style=style, mutate=mutate)
					self.assertAlmostEqual(PETRA.getExperimentInfo(path)["pix_size"], pix_size)
					images = PETRA.imageDataset(path, includeCurrent=False, includePixelShift=True)
					self.assertTrue(np.isfinite(images["pixel_shift"]).all())
		info = PETRA.getExperimentInfo(self.scan(mutate=old_names))
		self.assertAlmostEqual(info["field_of_view_x"], 130.0)  # whole sensor without a ROI
		self.assertEqual(info["camera"]["roi_width"], 200.0)

	def test_empty_strings_and_missing_groups(self):
		def empty_name(h5):
			replace(h5, "entry/beamline/name", np.array([], dtype=object), dtype=h5py.string_dtype())

		def no_setup(h5):
			del h5["entry/scan/setup"]
		info = PETRA.getExperimentInfo(self.scan(mutate=empty_name))
		self.assertNotIn("beamline", info)
		info = PETRA.getExperimentInfo(self.scan(mutate=no_setup))
		self.assertNotIn("energy_keV", info)
		self.assertAlmostEqual(info["pix_size"], 0.65)

	def test_incomplete_logs_are_skipped(self):
		def string_log(h5):
			data = h5["entry/scan/data"]
			data.create_dataset("comment/time", data=data["image_key/time"][:])
			data.create_dataset("comment/value", data=np.array([b"c"] * len(KEYS), dtype=object), dtype=h5py.string_dtype())

		def broken_logs(h5):
			data = h5["entry/scan/data"]
			data.create_dataset("no_value/time", data=data["image_key/time"][:3])
			data.create_dataset("not_a_group", data=np.arange(3))
			data.create_dataset("short_value/time", data=data["image_key/time"][:])
			data.create_dataset("short_value/value", data=np.arange(3.0))
		df = PETRA.scanDataset(self.scan(mutate=string_log))
		self.assertEqual(list(df["comment"]), ["c"] * len(KEYS))
		df = PETRA.scanDataset(self.scan(mutate=broken_logs))
		self.assertEqual(list(df.columns), ["time", "image_key", "image_file", "s_rot", "s_stage_x"])

	def test_truncated_image_file(self):
		def truncated(h5):  # an interrupted write: image_file shorter than image_key
			for part in ("value", "time"):
				v = h5[f"entry/scan/data/image_file/{part}"][:-3]
				replace(h5, f"entry/scan/data/image_file/{part}", v,
						dtype=h5py.string_dtype() if part == "value" else None)
		path = self.scan(mutate=truncated)
		info = PETRA.getExperimentInfo(path)  # from the image_key log
		self.assertEqual(info["sample_frame_count"], 4)
		with self.assertRaisesRegex(ValueError, "Length mismatch"):
			PETRA.scanDataset(path)

	def test_non_increasing_times(self):
		def duplicate(h5):
			for name in ("image_key", "image_file"):
				t = h5[f"entry/scan/data/{name}/time"][:]
				t[6] = t[5]
				replace(h5, f"entry/scan/data/{name}/time", t)
			del h5["entry/scan/data/s_rot"]
			del h5["entry/scan/data/s_stage_x"]
		df = PETRA.scanDataset(self.scan(mutate=duplicate), includeCurrent=True)
		self.assertTrue(np.all(np.diff(df.index.to_numpy().astype(np.int64)) > 0))
		self.assertTrue(isinstance(df["time"].iloc[0], pd.Timestamp))


class TestOneReader(unittest.TestCase):
	"""All functions read through one opening of the file and take an open h5py object; group
	attributes are values as the scalar datasets (as lautools PETRA, the canonical description)."""

	def test_open_file_and_attributes(self):
		def mutate(h5):
			h5["entry/hardware/camera"].attrs["model"] = np.array(["test camera"], dtype=object)
			h5["entry/hardware/camera"].attrs["sensor_gain"] = np.array([2.0])
			calibration = h5.create_group("entry/hardware/camera/calibration")
			calibration.attrs["magnification"] = 5.0
			del h5["entry/hardware/camera/magnification"]
			del h5["entry/beamline/experiment"]
			h5["entry/beamline"].attrs["experiment"] = "microtomography"
		with tempfile.TemporaryDirectory() as tmp:
			path = make_scan(os.path.join(tmp, "scan.h5"), style="2021", mutate=mutate)
			info = PETRA.getExperimentInfo(path)
			self.assertEqual(info["camera"]["model"], "test camera")
			self.assertEqual(info["camera"]["sensor_gain"], 2.0)
			self.assertEqual(info["camera"]["magnification"], 5.0)  # from the attribute of calibration
			self.assertEqual(info["experiment"], "microtomography")  # an attribute of beamline
			self.assertNotIn("NX_class", info["camera"])
			with h5py.File(path, "r") as h5:
				self.assertEqual(PETRA.getExperimentInfo(h5), info)
				pd.testing.assert_frame_equal(PETRA.scanDataset(h5, includeCurrent=True), PETRA.scanDataset(path, includeCurrent=True))
				pd.testing.assert_frame_equal(PETRA.imageDataset(h5["entry"], includePixelShift=True), PETRA.imageDataset(path, includePixelShift=True))
				pd.testing.assert_frame_equal(PETRA.getBasicDatasetInfo(h5), PETRA.getBasicDatasetInfo(path))
				self.assertTrue(bool(h5))  # not closed by the functions


if __name__ == "__main__":
	unittest.main()
