#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Processing of h5 format used to output Petra III data

@author: Vojtech Kulvait
@date: 2022-2025
"""
import h5py
import logging
import os
import pandas as pd
import numpy as np
import sys
from bisect import bisect_left
from scipy.signal import find_peaks

# Create a logger specific to this module
log = logging.getLogger(__name__)
log.setLevel(logging.INFO)  # Set the logging level to INFO
# Create a console handler and set its level to INFO
ch = logging.StreamHandler()
ch.setLevel(logging.INFO)
# Create a formatter and set it for the handler
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
ch.setFormatter(formatter)
# Add the handler to the logger
log.addHandler(ch)


def findInsertionEvents(h5file, timeOffsetSec=None, zeroTime=None):
	"""
	Finds insertion events in a PETRA III HDF5 file.

	Parameters:
			h5file (str or h5py.File): Path to the HDF5 file or an open HDF5 file object.
			timeOffsetSec (float, optional): Time offset in seconds to adjust timestamps from h5 file.
			zeroTime (int, optional): If specified, then start_time, end_time and mid_time will be returned as second offset relative to this value.

	Returns:
			List of insertion events with timestamps and values.
	"""
	if timeOffsetSec is None:
		timestampadjustment = np.int64(0)
	else:
		timestampadjustment = np.int64(timeOffsetSec * 1000)
	# Determine if h5file is a file path or already an h5py.File object
	h5_opened_here = False
	if isinstance(h5file, str):
		h5 = h5py.File(h5file, 'r')
		h5_opened_here = True
	else:
		h5 = h5file
	# Get beam current data with provided time offset
	beamCurrentData = beamCurrentDataset(h5, timeOffsetSec=timeOffsetSec)
	if zeroTime is not None:
		beamCurrentTimes = [(x - zeroTime).to_pytimedelta().total_seconds() for x in beamCurrentData["time"]]
	else:
		beamCurrentTimes = beamCurrentData["time"]
	beamCurrentValues = beamCurrentData["current"].to_numpy()
	# 1. Find local minima and maxima
	min_indices, _ = find_peaks(-beamCurrentValues, distance=5)
	max_indices, _ = find_peaks(beamCurrentValues, distance=5)
	# 2. Pair each min with next max
	insertion_events = []
	i = j = 0
	while i < len(min_indices) and j < len(max_indices):
		min_idx = min_indices[i]
		max_idx = max_indices[j]
		if max_idx > min_idx:
			# 3. Validate monotonic increase between min and max
			segment = beamCurrentValues[min_idx:max_idx + 1]
			if np.all(np.diff(segment) >= 0):  # monotonic non-decreasing
				# Optionally: filter by time or value thresholds
				if zeroTime is not None:
					dt = beamCurrentTimes[max_idx] - beamCurrentTimes[min_idx]
				else:
					dt = (beamCurrentTimes.iloc[max_idx] -
					      beamCurrentTimes.iloc[min_idx]).to_pytimedelta().total_seconds()
				dI = beamCurrentValues[max_idx] - beamCurrentValues[min_idx]
				if 1.0 < dt < 10.0 and dI > 0.5:  # tune these as needed
					insertion_events.append({
					    "start_index": min_idx,
					    "end_index": max_idx,
					    "start_time": beamCurrentTimes[min_idx],
					    "end_time": beamCurrentTimes[max_idx],
					    "mid_time": 0.5 * (beamCurrentTimes[min_idx] + beamCurrentTimes[max_idx]),
					    "current_increase": dI,
					    "duration_sec": dt
					})
			i += 1
		else:
			j += 1
	# Close the HDF5 file if it was opened here
	if h5_opened_here:
		h5.close()
	return insertion_events


# This function extracts beam current readings from PETRA III scan HDF5 files.
# Creates a DataFrame with time and current values in mA.
# These readings are used to relate frame time data with synchrotron beam current.
# The index is based on UNIX timestamps in milliseconds, allowing accurate interpolation.
# `timeOffsetSec` is applied to synchronize measurements with external frame time references.
def beamCurrentDataset(h5file, timeOffsetSec=None):
	if timeOffsetSec is None:
		timestampadjustment = np.int64(0)
	else:
		timestampadjustment = np.int64(timeOffsetSec * 1000)
	# Determine if h5file is a file path or already an h5py.File object
	h5_opened_here = False
	if isinstance(h5file, str):
		h5 = h5py.File(h5file, 'r')
		h5_opened_here = True
	else:
		h5 = h5file
	# Extract raw timestamp IDs (in milliseconds since epoch)
	ID = list(h5["entry/hardware/beam_current/current/time"])
	# Convert raw timestamps to datetime objects for human-readable interpretation
	time = list(pd.to_datetime(i, unit="ms") for i in ID)
	# Extract beam current values (in mA) from the HDF5 file
	value = list(h5["entry/hardware/beam_current/current/value"])
	# Combine time and current values into a list of tuples
	agg = list(zip(time, value))
	# Create a DataFrame from the combined data, with timestamps as the index
	df = pd.DataFrame(agg, columns=["time", "current"], index=ID)
	# Remove duplicate timestamps (keep first occurrence only)
	# This is important to avoid issues in interpolation or reindexing
	# see https://stackoverflow.com/questions/13035764/remove-pandas-rows-with-duplicate-indices
	df = df[~df.index.duplicated(keep='first')]
	# Remove any rows where the index (timestamp) is exactly zero — likely invalid
	# Remove zero index https://stackoverflow.com/questions/13851535/how-to-delete-rows-from-a-pandas-dataframe-based-on-a-conditional-expression
	df.drop(df[df.index == 0].index, inplace=True)
	# Sort by index (timestamp) to ensure chronological order
	df.sort_index(inplace=True)
	# Apply a time offset to the index (timestamp) to synchronize with external data
	if timestampadjustment != 0:
		# Adjustment in milliseconds
		df.index = (df.index + timestampadjustment).astype(np.uint64)
		# Update the time column to reflect the adjusted timestamps
		df.time = pd.to_datetime(df.index, unit="ms")
	# Close the HDF5 file if it was opened here
	if h5_opened_here:
		h5.close()
	return df


def getExperimentInfo(h5file, overrideMagnification=None):
	h5 = h5py.File(h5file, 'r')
	info = {}
	info["h5"] = os.path.abspath(h5file)
	setup = {}
	camera = {}
	if 'entry/scan/setup' in h5:
		setup_group = h5['entry/scan/setup']
		# Iterate over all items in the setup group
		for key in setup_group.keys():
			# Assuming each entry is a dataset containing a single double value
			setup[key] = setup_group[key][()][0]
	else:
		log.warning("The path 'entry/scan/setup' does not exist in %s." % info["h5"])
	# Check for the camera data path
	if 'entry/hardware/camera' in h5:
		camera_group = h5['entry/hardware/camera']
	elif 'entry/hardware/camera1' in h5:
		camera_group = h5['entry/hardware/camera1']
	else:
		camera_group = None
		log.warning("Neither 'entry/hardware/camera' nor 'entry/hardware/camera1' exists in %s." % info["h5"])
	# If a camera group was found, extract the data
	if camera_group is not None:
		for key in camera_group.keys():
			# Assuming each entry is a dataset containing a single double value
			camera[key] = camera_group[key][()][0]
		if overrideMagnification is not None:
			camera["magnification"] = overrideMagnification
		if "magnification" in camera and "pixelsize" in camera:
			info["pix_size"] = camera["pixelsize"] / camera["magnification"]
			info["pix_size_cam"] = camera["pixelsize"]
	info["setup"] = setup
	info["camera"] = camera
	return info


# Funcion no longer used, kept for reference use insert_label instead
def _insertToDf(df, dat, name):
	time = dat["%s/time" % (name)]
	value = dat["%s/value" % (name)]
	if name == "image_file" and h5py.__version__.startswith("3"):
		value = value.asstr()
	if len(time) != len(value):
		raise IOError("len(time)=%d len(value)=%d for entry %s" % (len(time), len(value), name))
	for i in range(len(value)):
		t = time[i]
		v = value[i]
		if t in df.index:
			df.loc[t][name] = v
		else:
			raise IOError("Can not insert value=%s into the column=%s because related time=%s is not in index" %
			              (v, name, t))


def normalize_timestamps_preserve_increase(old_ts):
	"""
	Normalize timestamps so that:
		new[i] = new[i-1] + max(1, old[i] - old[i-1])
	preserving original positive increases and fixing zero/negative increases.
	"""
	old_ts = np.asarray(old_ts, dtype=np.int64)
	if len(old_ts) == 0:
		return old_ts
	new_ts = old_ts.copy()
	for i in range(1, len(old_ts)):
		original_increase = old_ts[i] - old_ts[i - 1]
		required_increase = max(1, original_increase)
		new_ts[i] = new_ts[i - 1] + required_increase
	return new_ts


def getBasicDatasetInfo(h5file):
	with h5py.File(h5file, 'r') as h5:
		data = h5["entry/scan/data"]
		key_time = data["image_key/time"]
		key_value = data["image_key/value"]
		file_time = data["image_file/time"]
		file_value = data["image_file/value"]
		if h5py.__version__.startswith("3"):
			file_value = file_value.asstr()[:]
		# Check that all four arrays have the same length
		lengths = [len(key_time), len(key_value), len(file_time), len(file_value)]
		if len(set(lengths)) != 1:
			raise ValueError(f"Length mismatch: key_time={len(key_time)}, key_value={len(key_value)}, "
			                 f"file_time={len(file_time)}, file_value={len(file_value)}")
		if not np.array_equal(key_time, file_time):
			raise ValueError("Time arrays for image_key and image_file do not match.")
		# Create a DataFrame
		df = pd.DataFrame({"time": key_time, "image_key": key_value, "image_file": file_value})
		# Check for duplicate times
		duplicates = df[df.duplicated(subset=['time'], keep=False)].copy()
		if not duplicates.empty:
			# Save original indices
			duplicates["orig_index"] = duplicates.index
			n_duplicated_times = duplicates['time'].nunique()
			abs_h5file = os.path.abspath(h5file)
			print("WARNING: There is %d duplicated entries in %d distinct times in file %s." %
			      (len(duplicates), n_duplicated_times, abs_h5file))
			# Build key sets per time
			times_to_duplicates = duplicates.groupby('time')['image_key'].apply(list).to_dict()
			duplicates = duplicates.sort_values(by=['time', 'image_key'])
			duplicates["key_set"] = duplicates["time"].map(lambda t: tuple(sorted(set(times_to_duplicates[t]))))
			# Report per key_set
			for key_set, group in duplicates.groupby('key_set'):
				first_time = group['time'].iloc[0]
				affected = group[group['time'] == first_time][['image_file', 'orig_index']]
				formatted_key_set = ", ".join(str(k) for k in key_set)
				# Format as "file (id=index)"
				affected_str_list = [f"{f} (id={i})" for f, i in zip(affected['image_file'], affected['orig_index'])]
				affected_str = ", ".join(affected_str_list)
				print(f"key_set=[{formatted_key_set}], n_times={group['time'].nunique()}, "
				      f"first_time={first_time}, affected_files={affected_str}")
		return df


def insert_label(df, data, label_name, key_col="image_key"):
	"""
	Insert an HDF5 dataset 'label_name' into DataFrame 'df' based on timestamp alignment.
	Handles full DataFrame or specific filtered key subsets.
	"""
	time_ds = np.array(data[f"{label_name}/time"])
	value_ds = np.array(data[f"{label_name}/value"])
	if h5py.__version__.startswith("3") and value_ds.dtype.kind == 'O':
		value_ds = value_ds.asstr()[:]
	# Case 1: Full DataFrame
	if len(time_ds) == len(df):
		if not np.all(time_ds == df["time"].values):
			raise ValueError(f"{label_name}/time does not match df['time']")
		df[label_name] = value_ds
		return df
	# Case 2: Subsets based on key filters
	# Define candidate filters
	filters = [
	    (df[key_col] != 2),  # exclude key==2
	    ((df[key_col] != 2) & (df[key_col] != 1))  # exclude key 1 and 2
	]
	for mask in filters:
		if mask.sum() == len(time_ds):
			# Candidate subset length matches
			if not np.all(time_ds == df.loc[mask, "time"].values):
				raise ValueError(f"{label_name}/time shape matches subset ({mask.sum()} rows) "
				                 "but timestamps do not align exactly.")
			df.loc[mask, label_name] = value_ds
			return df
	# If no matching subset found, raise error
	raise ValueError(f"Could not insert {label_name}: no matching subset found.")


#Note that timeOffsetSec describes the offset to be applied to the current measurements to be in sync with the frame time data
def scanDataset(h5file, includeCurrent=False, timeOffsetSec=None):
	if timeOffsetSec is None:
		timeOffsetSec = 0.0
	h5 = h5py.File(h5file, 'r')
	abs_h5file = os.path.abspath(h5file)
	data = h5["entry/scan/data"]
	labels = list(data.keys())
	if len(labels) < 1:
		raise ValueError("Error: labels count is %d!" % (len(labels)))
	# Create basic dataset info based on image_key and image_file
	df = getBasicDatasetInfo(h5file)
	# Insert labels except image_key and image_file
	for lab in labels:
		if lab not in ["image_key", "image_file"]:
			try:
				df = insert_label(df, data, lab, key_col="image_key")
			except ValueError as ve:
				print(f"ValueError inserting label '{lab}' into DataFrame from file {abs_h5file}: {ve}")
				raise
	# Normalize timestamps to ensure strictly increasing index
	time_values = df["time"].values.astype(np.int64)
	if not all(np.diff(time_values) > 0):
		print("WARNING: Non-increasing timestamps detected in file %s; normalizing timestamps." % abs_h5file)
		normalized_time = normalize_timestamps_preserve_increase(df["time"].values)
		df["time"] = normalized_time
	# Set index to time for further processing
	df.index = df["time"].values
	# Convert time column to datetime in-place
	df["time"] = pd.to_datetime(df["time"], unit="ms")
	if includeCurrent:
		df["current"] = np.nan
		currentFrame = beamCurrentDataset(h5file, timeOffsetSec=timeOffsetSec)
		for ind in df.index:
			posAfterEq = bisect_left(currentFrame.index, ind)
			if posAfterEq == len(currentFrame.index):
				df.loc[ind, "current"] = currentFrame.iloc[posAfterEq - 1]["current"]
			elif posAfterEq == 0 or currentFrame.index[posAfterEq] == ind:
				df.loc[ind, "current"] = currentFrame.iloc[posAfterEq]["current"]
			else:
				t0 = float(currentFrame.index[posAfterEq - 1])
				cur0 = float(currentFrame.iloc[posAfterEq - 1]["current"])
				t1 = float(currentFrame.index[posAfterEq])
				cur1 = float(currentFrame.iloc[posAfterEq]["current"])
				t = float(ind)
				cur = cur0 + ((t - t0) / (t1 - t0)) * (cur1 - cur0)
				if np.isnan(cur):
					print("Cur is NaN t0=%f cur0=%f t1=%f cur1=%f t=%f ind=%d posAfterEq=%d len(currentFrame)=%d!" %
					      (t0, cur0, t1, cur1, t, ind, posAfterEq, len(currentFrame)))
					raise ValueError("Interpolation error producing NaN beam curent.")
				df.loc[ind, "current"] = cur
	df = df.sort_values("time", ascending=True)
	return df


def imageDataset(h5file, image_key=0, includeCurrent=True, includePixelShift=False, overrideMagnification=None):
	df = scanDataset(h5file, includeCurrent)
	df = df.loc[df["image_key"] == image_key]
	df = df.assign(frame_ind=np.arange(len(df)))
	if includePixelShift:
		h5 = h5py.File(h5file, 'r')
		if "/entry/hardware/camera1" in h5:
			cam = "camera1"
		elif "/entry/hardware/camera" in h5:
			cam = "camera"
		else:
			raise ValueError("There is no entry/hardware/camera or entry/hardware/camera1 entry in %s." % info["h5"])
		pix_size_cam = float(h5["entry/hardware/%s/pixelsize" % cam][0])
		if overrideMagnification is not None:
			pix_size_mag = overrideMagnification
		else:
			pix_size_mag = float(h5["entry/hardware/%s/magnification" % cam][0])
		pix_size = float(pix_size_cam / pix_size_mag)
		# zeropos = df["s_stage_x"].iloc[0]#Previous
		zeropos = 0.5 * (df["s_stage_x"].max() + df["s_stage_x"].min())
		# print("p0 =%f min=%f max=%f zeropos=%f pix_size=%f"%(df["s_stage_x"].iloc[0], df["s_stage_x"].min(), df["s_stage_x"].max(), zeropos, pix_size))
		pixShifts = (df["s_stage_x"] - zeropos) / pix_size
		df = df.assign(pixel_shift=pixShifts)
	return df
