#!/usr/bin/env python
#-*- coding: utf-8 -*-
"""
Processing of h5 format used to output Petra III data

@author: Vojtech Kulvait
@date: 2022-2024
"""
import h5py
import logging
import os
import pandas as pd
import numpy as np
from bisect import bisect_left

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


def beamCurrentDataset(h5file):
	h5 = h5py.File(h5file, 'r')
	ID = list(h5["entry/hardware/beam_current/current/time"])
	time = list(pd.to_datetime(i, unit="ms") for i in ID)
	value = list(h5["entry/hardware/beam_current/current/value"])
	agg = list(zip(time, value))
	df = pd.DataFrame(agg, columns=["time", "current"], index=ID)
	#There were duplicates caused problems
	#https://stackoverflow.com/questions/13035764/remove-pandas-rows-with-duplicate-indices
	df = df[~df.index.duplicated(keep='first')]
	#Remove zero index https://stackoverflow.com/questions/13851535/how-to-delete-rows-from-a-pandas-dataframe-based-on-a-conditional-expression
	df.drop(df[df.index == 0].index, inplace=True)
	#Need to sort it for interpolation to work
	df.sort_index(inplace=True)
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


def scanDataset(h5file, includeCurrent=False):
	h5 = h5py.File(h5file, 'r')
	data = h5["entry/scan/data"]
	labels = list(data.keys())
	if len(labels) < 1:
		sys.exit("Error: labels count is %d!" % (labels.count))
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
		currentFrame = beamCurrentDataset(h5file)
		for ind in df.index:
			posAfterEq = bisect_left(currentFrame.index, ind)
			if posAfterEq == len(currentFrame.index):
				df.loc[ind]["current"] = currentFrame.iloc[posAfterEq -
														   1]["current"]
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
