#!/usr/bin/env python
#-*- coding: utf-8 -*-
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
log.setLevel(logging.INFO)	# Set the logging level to INFO
# Create a console handler and set its level to INFO
ch = logging.StreamHandler()
ch.setLevel(logging.INFO)
# Create a formatter and set it for the handler
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
ch.setFormatter(formatter)
# Add the handler to the logger
log.addHandler(ch)

def _insertToDf(df, dat, name):
	time = dat["%s/time" % (name)]
	value = dat["%s/value" % (name)]
	if name == "image_file" and h5py.__version__.startswith("3"):
		value = value.asstr()
	if len(time) != len(value):
		raise IOError("len(time)=%d len(value)=%d for entry %s" %
					  (len(time), len(value), name))
	for i in range(len(value)):
		t = time[i]
		v = value[i]
		if t in df.index:
			df.loc[t][name] = v
		else:
			raise IOError(
				"Can not insert value=%s into the column=%s because related time=%s is not in index"
				% (v, name, t))


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
			segment = beamCurrentValues[min_idx:max_idx+1]
			if np.all(np.diff(segment) >= 0):  # monotonic non-decreasing
				# Optionally: filter by time or value thresholds
				if zeroTime is not None:
					dt = beamCurrentTimes[max_idx] - beamCurrentTimes[min_idx]
				else:
					dt = (beamCurrentTimes.iloc[max_idx] - beamCurrentTimes.iloc[min_idx]).to_pytimedelta().total_seconds()
				dI = beamCurrentValues[max_idx] - beamCurrentValues[min_idx]
				if 1.0 < dt < 60.0 and dI > 0.05:  # tune these as needed
					insertion_events.append({
						"start_index": min_idx,
						"end_index": max_idx,
						"start_time": beamCurrentTimes[min_idx],
						"end_time": beamCurrentTimes[max_idx],
						"mid_time": 0.5*(beamCurrentTimes[min_idx]+beamCurrentTimes[max_idx]),
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
	#Combine time and current values into a list of tuples
	agg = list(zip(time, value))
	# Create a DataFrame from the combined data, with timestamps as the index
	df = pd.DataFrame(agg, columns=["time", "current"], index=ID)
	# Remove duplicate timestamps (keep first occurrence only)
	# This is important to avoid issues in interpolation or reindexing
	# see https://stackoverflow.com/questions/13035764/remove-pandas-rows-with-duplicate-indices
	df = df[~df.index.duplicated(keep='first')]
	# Remove any rows where the index (timestamp) is exactly zero — likely invalid
	#Remove zero index https://stackoverflow.com/questions/13851535/how-to-delete-rows-from-a-pandas-dataframe-based-on-a-conditional-expression
	df.drop(df[df.index == 0].index, inplace=True)
	# Sort by index (timestamp) to ensure chronological order
	df.sort_index(inplace=True)
	# Apply a time offset to the index (timestamp) to synchronize with external data
	if timestampadjustment != 0:
		df.index = (df.index + timestampadjustment).astype(np.uint64) # Adjustment in milliseconds
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

#Note that timeOffsetSec describes the offset to be applied to the current measurements to be in sync with the frame time data
def scanDataset(h5file, includeCurrent=False, timeOffsetSec=None):
	if timeOffsetSec is None:
		timeOffsetSec = 0.0
	h5 = h5py.File(h5file, 'r')
	data = h5["entry/scan/data"]
	labels = list(data.keys())
	if len(labels) < 1:
		raise ValueError("Error: labels count is %d!" % (len(labels)))
	#There always shall be image_key entry
	if includeCurrent:
		df = pd.DataFrame(columns=labels + ["current", "time"],
						  index=list(data["image_key/time"]))
	else:
		df = pd.DataFrame(columns=labels + ["time"],
						  index=list(data["image_key/time"]))
	for ind in df.index:
		df.loc[ind, "time"] = pd.to_datetime(ind, unit="ms")
	for lab in labels:
		try:
			_insertToDf(df, data, lab)
		except IOError:
			print("IOError processing %s" % h5file)
			raise
	if includeCurrent:
		currentFrame = beamCurrentDataset(h5file, timeOffsetSec=timeOffsetSec)
		for ind in df.index:
			posAfterEq = bisect_left(currentFrame.index, ind)
			if posAfterEq == len(currentFrame.index):
				df.loc[ind]["current"] = currentFrame.iloc[posAfterEq - 1]["current"]
			elif posAfterEq == 0 or currentFrame.index[posAfterEq] == ind:
				df.loc[ind]["current"] = currentFrame.iloc[posAfterEq][
					"current"]
			else:
				t0 = float(currentFrame.index[posAfterEq - 1])
				cur0 = float(currentFrame.iloc[posAfterEq - 1]["current"])
				t1 = float(currentFrame.index[posAfterEq])
				cur1 = float(currentFrame.iloc[posAfterEq]["current"])
				t = float(ind)
				cur = cur0 + ((t - t0) / (t1 - t0)) * (cur1 - cur0)
				if np.isnan(cur):
					print(
						"Cur is NaN t0=%f cur0=%f t1=%f cur1=%f t=%f ind=%d posAfterEq=%d len(currentFrame)=%d!"
						% (t0, cur0, t1, cur1, t, ind, posAfterEq,
						   len(currentFrame)))
					raise ValueError(
						"Interpolation error producing NaN beam curent.")
				df.loc[ind]["current"] = cur
	df = df.sort_values("time", ascending=True)
	return df

def imageDataset(h5file,
				 image_key=0,
				 includeCurrent=True,
				 includePixelShift=False,
				 overrideMagnification=None):
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
			raise ValueError(
				"There is no entry/hardware/camera or entry/hardware/camera1 entry in %s."
				% info["h5"])
		pix_size_cam = float(h5["entry/hardware/%s/pixelsize" % cam][0])
		if overrideMagnification is not None:
			pix_size_mag = overrideMagnification
		else:
			pix_size_mag = float(h5["entry/hardware/%s/magnification" % cam][0])
		pix_size = float(pix_size_cam / pix_size_mag)
		#zeropos = df["s_stage_x"].iloc[0]#Previous
		zeropos = 0.5 * (df["s_stage_x"].max() + df["s_stage_x"].min())
		#print("p0 =%f min=%f max=%f zeropos=%f pix_size=%f"%(df["s_stage_x"].iloc[0], df["s_stage_x"].min(), df["s_stage_x"].max(), zeropos, pix_size))
		pixShifts = (df["s_stage_x"] - zeropos) / pix_size
		df = df.assign(pixel_shift=pixShifts)
	return df
