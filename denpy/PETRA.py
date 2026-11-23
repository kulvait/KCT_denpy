#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Processing of h5 format used to output Petra III data

@author: Vojtech Kulvait
@date: 2022-2026
"""
import contextlib
import h5py
import os
import pandas as pd
import numpy as np
import sys
from scipy.signal import find_peaks

import logging
# Create a logger specific to this module
log = logging.getLogger(__name__)
log.setLevel(logging.INFO) # Set the logging level to INFO
# Create a console handler and set its level to INFO
ch = logging.StreamHandler()
ch.setLevel(logging.INFO)
# Create a formatter and set it for the handler
formatter = logging.Formatter('%(asctime)s - %(name)s:%(lineno)d - %(levelname)s : %(message)s', datefmt='%d.%m.%Y %H:%M:%S')
ch.setFormatter(formatter)
# Add the handler to the logger
log.addHandler(ch)
log.propagate = False # Prevent log messages from being propagated to the root logger


class BeamCurrentError(ValueError):
	"""The beam current of a scan file is missing or unusable (entry/hardware/beam_current/current)."""


@contextlib.contextmanager
def _open_h5(h5file):
	"""The HDF5 file of a file name (opened and closed here), or the file of an open h5py.File or
	h5py.Group (left open). All functions of this module read through it, so a function opens the
	file once, also for the functions it calls."""
	if isinstance(h5file, (str, os.PathLike)):
		with h5py.File(h5file, "r") as h5:
			yield h5
	else:
		yield h5file.file


def _scalar_number(value):
	"""value as a finite float if it is a single number (not bool), else None."""
	if isinstance(value, np.ndarray) and value.size == 1:
		value = value.reshape(()).item()
	if isinstance(value, np.generic):
		value = value.item()
	if isinstance(value, bool) or not isinstance(value, (int, float)):
		return None
	value = float(value)
	return value if np.isfinite(value) else None


def findInsertionEvents(
	h5file,
	timeOffsetSec=None,
	zeroTime=None,
	*,
	peakDistance=5,
	minDurationSec=1.0,
	maxDurationSec=10.0,
	minCurrentIncrease=0.5,
):
	"""Find current insertion events from H5 metadata or a DataFrame.

	DataFrame input must contain datetime 'time' and numeric 'current'
	columns. Event indices refer to row positions in that DataFrame.

	With zeroTime, event times are seconds relative to zeroTime.
	Otherwise, event times are pandas timestamps.
	"""
	offset = 0.0 if timeOffsetSec is None else float(timeOffsetSec)
	if not np.isfinite(offset):
		raise ValueError("timeOffsetSec must be finite.")
	if peakDistance < 1:
		raise ValueError("peakDistance must be at least 1.")
	if not 0 <= minDurationSec < maxDurationSec:
		raise ValueError("Invalid insertion duration bounds.")
	if not np.isfinite(minCurrentIncrease) or minCurrentIncrease < 0:
		raise ValueError("minCurrentIncrease must be finite and nonnegative.")
	if isinstance(h5file, pd.DataFrame):
		missing = {"time", "current"} - set(h5file.columns)
		if missing:
			raise ValueError(f"Current dataset is missing columns: {missing}")
		data = h5file[["time", "current"]].copy()
		# Do not modify the caller's DataFrame.
		if pd.api.types.is_numeric_dtype(data["time"].dtype):
			raise TypeError("DataFrame time values must be datetimes.")
		data["time"] = pd.to_datetime(data["time"], utc=True, errors="raise")
		data["time"] += pd.to_timedelta(offset, unit="s")
	else:
		# beamCurrentDataset owns/cleans up any file it opens.
		source = os.fspath(h5file) if isinstance(h5file, os.PathLike) else h5file
		data = beamCurrentDataset(source, timeOffsetSec=offset).copy()
		data["time"] = pd.to_datetime( data["time"], utc=True, errors="raise")
	data = data.reset_index(drop=True)
	values = pd.to_numeric(
		data["current"], errors="raise"
	).to_numpy(dtype=np.float64)
	if data["time"].isna().any() or not np.all(np.isfinite(values)):
		raise ValueError("Current timestamps and values must be valid.")
	if len(data) < 2:
		return []
	if not (data["time"].diff().iloc[1:] > pd.Timedelta(0)).all():
		raise ValueError("Current timestamps must be strictly increasing.")
	times = list(data["time"])
	origin = (
		None if zeroTime is None
		else pd.to_datetime(zeroTime, utc=True, errors="raise")
	)
	if origin is not None and pd.isna(origin):
		raise ValueError("zeroTime must be a valid datetime.")
	min_indices, _ = find_peaks(-values, distance=peakDistance)
	max_indices, _ = find_peaks(values, distance=peakDistance)
	events = []
	i = j = 0
	while i < len(min_indices) and j < len(max_indices):
		start = int(min_indices[i])
		end = int(max_indices[j])
		if end <= start:
			j += 1
			continue
		duration = (times[end] - times[start]).total_seconds()
		increase = float(values[end] - values[start])
		segment = values[start:end + 1]
		if (
			np.all(np.diff(segment) >= 0)
			and minDurationSec < duration < maxDurationSec
			and increase > minCurrentIncrease
		):
			if origin is None:
				start_time = times[start]
				end_time = times[end]
			else:
				start_time = (times[start] - origin).total_seconds()
				end_time = (times[end] - origin).total_seconds()
			events.append({
				"start_index": start,
				"end_index": end,
				"start_time": start_time,
				"end_time": end_time,
				"mid_time": start_time + (end_time - start_time) / 2,
				"current_increase": increase,
				"duration_sec": duration,
			})
		i += 1
	return events


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
	with _open_h5(h5file) as h5:
		filename = os.path.realpath(h5.filename)
		beam_current_path = "entry/hardware/beam_current/current"
		for part in ("time", "value"):
			if f"{beam_current_path}/{part}" not in h5:
				raise BeamCurrentError(f"No beam current in {filename}: '{beam_current_path}/{part}' does not exist.")
		# Raw timestamps (milliseconds since epoch) and beam current values (mA)
		ID = np.asarray(h5[f"{beam_current_path}/time"][()]).astype(np.uint64).ravel()
		value = np.asarray(h5[f"{beam_current_path}/value"][()], dtype=np.float64).ravel()
		if len(ID) != len(value):
			raise BeamCurrentError(f"Beam current in {filename} has {len(ID)} timestamps but {len(value)} values.")
		# Readings without a valid value carry no information
		valid = np.isfinite(value)
		if not valid.all():
			log.warning("Ignoring %d of %d beam current readings without a finite value in %s." % ((~valid).sum(), len(value), filename))
			ID, value = ID[valid], value[valid]
		# DataFrame with timestamps as the index
		df = pd.DataFrame({"time": pd.to_datetime(ID.astype(np.int64), unit="ms"), "current": value}, index=ID)
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
			# Apply the offset using signed arithmetic, independently of pandas/NumPy mixed signed/unsigned behavior.
			timestamps = df.index.to_numpy(dtype=np.int64)
			adjusted = timestamps + timestampadjustment
			# Check if uint64 Unix-millisecond index can be preserved
			if np.any(adjusted < 0):
				raise ValueError(f"Time offset {timeOffsetSec}s seconds results in negative timestamps. Invalid uint64 index.")
			# Adjustment in milliseconds
			df.index = pd.Index(adjusted.astype(np.uint64), name=df.index.name)
			# Update the time column to reflect the adjusted timestamps
			df["time"] = pd.to_datetime(df.index, unit="ms")
		return df


def read_scalar_h5(item):
	"""
	Read a scalar value from an h5py dataset robustly.

	Handles:
	  - shape ()		   -> 0-d array / numpy scalar
	  - shape (1,)		   -> 1-element array
	  - bytes / np.bytes_  -> decoded to str (utf-8, replace on error)
	  - numpy scalars	   -> converted to native Python types
	Returns None for empty datasets.
	"""
	val = item[()]
	# Unwrap 1-element arrays; keep 0-d as-is
	if isinstance(val, np.ndarray):
		if val.shape == (1,):
			val = val[0]
		elif val.shape == ():
			val = val.item()
		elif val.size == 0:
			return None
		# else: leave as ndarray (caller may want the array)
	# numpy scalar -> Python scalar (np.bytes_.item() returns bytes, handled below)
	if isinstance(val, np.generic):
		val = val.item()
	# Decode bytes to str
	if isinstance(val, (bytes, bytearray)):
		try:
			val = val.decode("utf-8")
		except UnicodeDecodeError:
			val = val.decode("utf-8", errors="replace")
	return val

def _attribute_value(value):
	"""An HDF5 attribute as a python value: numpy to python, bytes to str, arrays to lists; one-element
	sequences collapse to their element (older files store attributes as (1,) object arrays)."""
	if isinstance(value, np.ndarray):
		return _attribute_value(value.tolist())
	if isinstance(value, np.generic):
		return _attribute_value(value.item())
	if isinstance(value, (bytes, bytearray)):
		try:
			return value.decode("utf-8")
		except UnicodeDecodeError:
			return value.decode("utf-8", errors="replace")
	if isinstance(value, (list, tuple)):
		if len(value) == 1:
			return _attribute_value(value[0])
		return [_attribute_value(v) for v in value]
	return value


def _attribute_of_parent(h5, path):
	"""(True, value) if path names an attribute of its parent group (other than NX_class), else (False, None)."""
	parent, _, name = path.rpartition("/")
	if name != "NX_class" and (parent == "" or isinstance(h5.get(parent), h5py.Group)):
		group = h5[parent] if parent else h5
		if name in group.attrs:
			return True, _attribute_value(group.attrs[name])
	return False, None


#Sometimes values as entry/beamline/experiment are stored as () shape but sometimes as (1,) shape,
#or as an attribute of the group; this function extracts the string value from any of these.
def extractStringEntry(h5, path):
	is_attribute, value = _attribute_of_parent(h5, path)
	if path in h5 or is_attribute:
		if not is_attribute:
			value = read_scalar_h5(h5[path])
		if value is None:
			log.warning("The path '%s' is empty in %s." % (path, os.path.realpath(h5.filename)))
			return None
		if isinstance(value, str):
			return value
		else:
			return str(value)
	else:
		log.warning("The path '%s' does not exist in %s." % (path, os.path.realpath(h5.filename)))
		return None

def read_group_scalars(group):
	"""The datasets (read_scalar_h5) and the attributes other than NX_class of a group as a dict; an
	attribute wins over a dataset of the same name. Subgroups are skipped."""
	out = {}
	for key, item in group.items():
		if isinstance(item, h5py.Dataset):
			out[key] = read_scalar_h5(item)
	out.update({key: _attribute_value(value) for key, value in group.attrs.items() if key != "NX_class"})
	return out

def getExperimentInfo(h5file, overrideMagnification=None):
	with _open_h5(h5file) as h5:
		return _experimentInfo(h5, h5.filename, overrideMagnification)


def _experimentInfo(h5, h5file, overrideMagnification=None):
	info = {}
	info["h5"] = os.path.realpath(h5file)
	info["h5_name"] = os.path.basename(info["h5"])
	experimentString = extractStringEntry(h5, 'entry/beamline/experiment')
	if experimentString is not None:
		info["experiment"] = experimentString
	beamlineString = extractStringEntry(h5, 'entry/beamline/name')
	if beamlineString is not None:
		info["beamline"] = beamlineString
	# Extract start and end times from the image_file dataset if available
	if 'entry/scan/data/image_file/time' in h5:
		image_file_time = pd.to_datetime(np.asarray(h5['entry/scan/data/image_file/time'][:]).astype(np.int64), unit="ms")
		if len(image_file_time) > 0:
			start = image_file_time.min()
			end = image_file_time.max()
			info["start_date_time"] = start.strftime("%d.%m.%Y %H:%M:%S")
			info["end_date_time"] = end.strftime("%d.%m.%Y %H:%M:%S")
			info["scan_frame_count"] = len(image_file_time)
			info["scan_duration_sec"] = (end - start).total_seconds()
			info["scan_fps"] = len(image_file_time) / info["scan_duration_sec"] if info["scan_duration_sec"] > 0 else None
	# Sample frames (image_key 0) from the image_key log, its time and value have the same length
	if 'entry/scan/data/image_key/value' in h5 and 'entry/scan/data/image_key/time' in h5:
		image_key = np.asarray(h5['entry/scan/data/image_key/value'][:])
		image_key_time = np.asarray(h5['entry/scan/data/image_key/time'][:])
		if len(image_key) != len(image_key_time):
			log.warning("image_key has %d values but %d times in %s; sample frame statistics skipped." % (len(image_key), len(image_key_time), info["h5"]))
		else:
			image_times = pd.to_datetime(image_key_time[image_key == 0].astype(np.int64), unit="ms")
			if len(image_times) > 0:
				info["sample_frame_count"] = len(image_times)
				info["sample_duration_sec"] = (image_times.max() - image_times.min()).total_seconds()
				info["sample_fps"] = len(image_times) / info["sample_duration_sec"] if info["sample_duration_sec"] > 0 else None
	setup = {}
	camera = {}
	if 'entry/scan/setup' in h5:
		# Iterate over all items in the setup group
		setup = read_group_scalars(h5['entry/scan/setup'])
		for key in ("pos_o_ccd_dist", "o_ccd_dist"):
			distance = _scalar_number(setup.get(key))
			if distance is not None:
				info["propagation_distance_mm"] = round(distance, 0)
				break
		for key in ("pos_p05_energy", "pos_p07_energy", "p05_energy", "p07_energy"):
			energy = _scalar_number(setup.get(key))
			if energy is not None:
				info["energy_keV"] = round(energy / 1000, 3)
				break
	else:
		log.warning("The path 'entry/scan/setup' does not exist in %s." % info["h5"])
	camera_group = _camera_group(h5)
	if camera_group is None:
		log.warning("Neither 'entry/hardware/camera' nor 'entry/hardware/camera1' exists in %s." % info["h5"])
	else:
		camera = read_group_scalars(camera_group)
		values = _camera_values(camera_group, camera, overrideMagnification)
		if overrideMagnification is not None:
			camera["magnification"] = overrideMagnification
		elif values["magnification"] is not None:
			camera["magnification"] = values["magnification"]
		# Names of older files, as the newer ones
		for key in ("pixelsize", "sensorsize_x", "sensorsize_y", "roi_width", "roi_height"):
			if values[key] is not None:
				camera[key] = values[key]
		# Without a ROI, the whole sensor
		for roi, sensor in (("roi_width", "sensorsize_x"), ("roi_height", "sensorsize_y")):
			if values[roi] is None and values[sensor] is not None:
				camera[roi] = values[sensor]
		if values["magnification"] is not None and values["pixelsize"] is not None:
			info["pix_size"] = values["pixelsize"] / values["magnification"]
			width = values["roi_width"] if values["roi_width"] is not None else values["sensorsize_x"]
			height = values["roi_height"] if values["roi_height"] is not None else values["sensorsize_y"]
			if width is not None:
				info["field_of_view_x"] = info["pix_size"] * width
			if height is not None:
				info["field_of_view_y"] = info["pix_size"] * height
		else:
			log.warning("No usable camera pixelsize and magnification in %s; pix_size not set." % info["h5"])
	# Compute fresnel number if propagation distance, energy, and pixel size are available
	if "propagation_distance_mm" in info and "energy_keV" in info and "pix_size" in info and  info["energy_keV"] > 0 and info["propagation_distance_mm"] > 0 and info["pix_size"] > 0:
		wavelength_m = 1.2398419843320026e-9 / info["energy_keV"]  # Convert keV to eV for wavelength calculation
		propagation_distance_m = info["propagation_distance_mm"] / 1000.0  # Convert mm to m
		pix_size_m = info["pix_size"] / 1000.0	# Convert mm to m
		fresnel_number = (pix_size_m ** 2) / (wavelength_m * propagation_distance_m)
		info["fresnel_number"] = fresnel_number
	info["setup"] = setup
	info["camera"] = camera
	return info


def _camera_group(h5):
	"""entry/hardware/camera, else entry/hardware/camera1, else None."""
	for path in ("entry/hardware/camera", "entry/hardware/camera1"):
		if path in h5:
			return h5[path]
	return None


def _camera_values(camera_group, camera, overrideMagnification=None):
	"""Camera numbers under the names of the newer files, None where missing, empty or not a number.

	Older names: px_size/pixel_size, senzor_xsize/senzor_ysize, roi_xsize/roi_ysize. The
	magnification is taken from overrideMagnification, the camera group or its calibration subgroup;
	a magnification of 0 counts as missing.
	"""
	def first(*keys):
		for key in keys:
			value = _scalar_number(camera.get(key))
			if value is not None:
				return value
		return None
	values = {
		"pixelsize": first("pixel_size", "px_size", "pixelsize"),
		"sensorsize_x": first("senzor_xsize", "sensorsize_x"),
		"sensorsize_y": first("senzor_ysize", "sensorsize_y"),
		"roi_width": first("roi_xsize", "roi_width"),
		"roi_height": first("roi_ysize", "roi_height"),
	}
	magnification = _scalar_number(overrideMagnification) if overrideMagnification is not None else first("magnification")
	if magnification is None and overrideMagnification is None and isinstance(camera_group.get("calibration"), h5py.Group):
		magnification = _scalar_number(read_group_scalars(camera_group["calibration"]).get("magnification"))
	values["magnification"] = magnification if magnification else None
	return values


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
	with _open_h5(h5file) as h5:
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
			abs_h5file = os.path.realpath(h5.filename)
			log.warning("There is %d duplicated entries in %d distinct times in file %s." %
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
				log.warning(f"key_set=[{formatted_key_set}], n_times={group['time'].nunique()}, "
					  f"first_time={first_time}, affected_files={affected_str}")
		return df


def insert_label(df, data, label_name, key_col="image_key"):
	"""
	Insert an HDF5 dataset 'label_name' into DataFrame 'df' based on timestamp alignment.
	Handles full DataFrame or specific filtered key subsets.
	"""
	for part in ("time", "value"):
		if not isinstance(data.get(f"{label_name}/{part}"), h5py.Dataset):
			raise ValueError(f"{label_name} has no {part} dataset.")
	time_ds = np.asarray(data[f"{label_name}/time"][()])
	value_dataset = data[f"{label_name}/value"]
	if h5py.check_string_dtype(value_dataset.dtype) is not None:
		value_ds = np.asarray(value_dataset.asstr()[()], dtype=object)
	else:
		value_ds = np.asarray(value_dataset[()])
	if time_ds.ndim != 1 or value_ds.ndim != 1:
		raise ValueError(f"{label_name} time and value must be 1-D, got shapes {time_ds.shape} and {value_ds.shape}.")
	if len(time_ds) != len(value_ds):
		raise ValueError(f"{label_name} has {len(time_ds)} times but {len(value_ds)} values.")
	# In the situation time_ds and value_ds are empty, raise ValueError
	if len(time_ds) == 0:
		raise ValueError(f"{label_name} has empty time and value datasets.")
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
	"""DataFrame of the frames of a scan: time, image_key, image_file and the logs of entry/scan/data.

	With includeCurrent, the column current holds the beam current interpolated at the frame times
	(constant before the first and after the last reading); raises BeamCurrentError if the file has
	no usable beam current.
	"""
	if timeOffsetSec is None:
		timeOffsetSec = 0.0
	with _open_h5(h5file) as h5:
		return _scanDataset(h5, includeCurrent, timeOffsetSec)


def _scanDataset(h5, includeCurrent, timeOffsetSec):
	abs_h5file = os.path.realpath(h5.filename)
	# Create basic dataset info based on image_key and image_file
	df = getBasicDatasetInfo(h5)
	data = h5["entry/scan/data"]
	labels = list(data.keys())
	if len(labels) < 1:
		raise ValueError("Error: labels count is %d!" % (len(labels)))
	# Insert labels except image_key and image_file
	for lab in labels:
		if lab not in ["image_key", "image_file"]:
			if not isinstance(data[lab], h5py.Group):
				log.warning(f"Label '{lab}' not inserted into DataFrame from file {abs_h5file}: not a group with time and value.")
				continue
			try:
				df = insert_label(df.copy(deep=True), data, lab, key_col="image_key")  # Test on a copy to avoid partial insertions
			except ValueError as ve:
				log.warning(f"Label '{lab}' not inserted into DataFrame from file {abs_h5file}: {ve}")
	# Normalize timestamps to ensure strictly increasing index
	time_values = df["time"].values.astype(np.int64)
	if not all(np.diff(time_values) > 0):
		log.warning("Non-increasing timestamps detected in file %s; normalizing timestamps." % abs_h5file)
		normalized_time = normalize_timestamps_preserve_increase(df["time"].values)
		df["time"] = normalized_time
	# Set index to time for further processing
	df.index = df["time"].values
	# Convert time column to datetime in-place
	df["time"] = pd.to_datetime(df["time"].astype(np.int64), unit="ms")
	if includeCurrent:
		currentFrame = beamCurrentDataset(h5, timeOffsetSec=timeOffsetSec)
		if currentFrame.empty:
			raise BeamCurrentError(f"No beam current readings in {abs_h5file}, can not include the current.")
		current_times = currentFrame.index.to_numpy().astype(np.float64)
		frame_times = df.index.to_numpy().astype(np.float64)
		outside = int(np.sum(frame_times < current_times[0]) + np.sum(frame_times > current_times[-1]))
		if outside > 0:
			log.warning("%d of %d frames of %s are outside the beam current readings; the first/last reading is used for them." % (outside, len(df), abs_h5file))
		# Linear interpolation between the readings, constant outside of them
		df["current"] = np.interp(frame_times, current_times, currentFrame["current"].to_numpy(dtype=np.float64))
	df = df.sort_values("time", ascending=True)
	return df


def imageDataset(h5file, image_key=0, includeCurrent=True, includePixelShift=False, overrideMagnification=None):
	with _open_h5(h5file) as h5:
		return _imageDataset(h5, image_key, includeCurrent, includePixelShift, overrideMagnification)


def _imageDataset(h5, image_key, includeCurrent, includePixelShift, overrideMagnification):
	h5file = os.path.realpath(h5.filename)
	df = scanDataset(h5, includeCurrent)
	df = df.loc[df["image_key"] == image_key]
	df = df.assign(frame_ind=np.arange(len(df)))
	if includePixelShift:
		if "s_stage_x" not in df.columns:
			raise ValueError("No s_stage_x log in %s, can not compute the pixel shift." % h5file)
		camera_group = _camera_group(h5)
		if camera_group is None:
			raise ValueError("There is no entry/hardware/camera or entry/hardware/camera1 entry in %s." % h5file)
		values = _camera_values(camera_group, read_group_scalars(camera_group), overrideMagnification)
		if values["pixelsize"] is None or values["magnification"] is None:
			raise ValueError("No usable camera pixelsize and magnification in %s, can not compute the pixel shift." % h5file)
		pix_size = values["pixelsize"] / values["magnification"]
		# zeropos = df["s_stage_x"].iloc[0]#Previous
		zeropos = 0.5 * (df["s_stage_x"].max() + df["s_stage_x"].min())
		# print("p0 =%f min=%f max=%f zeropos=%f pix_size=%f"%(df["s_stage_x"].iloc[0], df["s_stage_x"].min(), df["s_stage_x"].max(), zeropos, pix_size))
		pixShifts = (df["s_stage_x"] - zeropos) / pix_size
		df = df.assign(pixel_shift=pixShifts)
	return df

def scanDatasetNanotomography(logfile, folder, includeCurrent=False):
	"""
	Build a scan DataFrame from LogScan.log and a directory of TIFF files.

	Files must follow:
		img_<batch>_<frame>.tif[f]
		ref_<batch>_<frame>.tif[f]

	Within each identifier, files are sorted numerically by (batch, frame)
	and paired with log rows in their original order. Log image numbers
	are preserved as metadata but are not used to construct filenames.

	Only log rows whose infostr is 'start' are used.

	Returns:
		DataFrame with an increasing UNIX-millisecond index and columns:
		time, image_key, image_file, s_rot, time_sec,
		image_identifier, infostr, image_number_i, image_number_ii,
		current_number, and optionally current.

	image_key:
		img -> 0 (projection)
		ref -> 1 (reference/flat field)

	No timeOffsetSec is needed: beam current is recorded in each log row.
	"""
	import io
	import re
	from decimal import Decimal

	folder = os.path.abspath(os.fspath(folder))
	logfile = os.path.abspath(os.fspath(logfile))
	if not os.path.isdir(folder):
		raise NotADirectoryError(folder)

	start_time = None
	rows = []
	with open(logfile, "r", encoding="utf-8") as stream:
		for line_number, line in enumerate(stream, 1):
			line = line.strip()
			if not line:
				continue
			if line.startswith("#"):
				match = re.fullmatch(
					r"#\s*starttime\s*=\s*(\S+)", line
				)
				if match:
					if start_time is not None:
						raise ValueError("Multiple #starttime headers.")
					start_time = Decimal(match.group(1))
				continue
			if len(line.split()) != 8:
				raise ValueError(
					f"{logfile}:{line_number}: expected 8 columns."
				)
			rows.append(line)

	if start_time is None or not start_time.is_finite():
		raise ValueError(f"Missing or invalid #starttime in {logfile}.")
	if not rows:
		raise ValueError(f"No data rows in {logfile}.")

	columns = [
		"image_identifier", "infostr",
		"image_number_i", "image_number_ii", "current_number",
		"time_sec", "current", "s_rot",
	]
	df = pd.read_csv(
		io.StringIO("\n".join(rows)),
		sep=r"\s+",
		names=columns,
		dtype=str,
	)
	df = df.loc[df["infostr"] == "start"].copy().reset_index(drop=True)
	if df.empty:
		raise ValueError(f"No 'start' rows in {logfile}.")

	key_map = {"img": 0, "ref": 1}
	unknown = set(df["image_identifier"]) - set(key_map)
	if unknown:
		raise ValueError(f"Unsupported image identifiers: {sorted(unknown)}.")

	# Convert epoch seconds + relative seconds directly to milliseconds.
	# Decimal avoids first rounding the epoch through a Python float.
	relative_times = [Decimal(value) for value in df["time_sec"]]
	if any(not value.is_finite() for value in relative_times):
		raise ValueError("Log timestamps must be finite.")
	timestamps = np.asarray(
		[int((start_time + value) * 1000) for value in relative_times],
		dtype=np.int64,
	)

	for column in ("image_number_i", "image_number_ii", "current_number"):
		df[column] = pd.to_numeric(df[column], errors="raise")
		if not np.isfinite(df[column].to_numpy()).all():
			raise ValueError(f"Non-finite values in {column}.")
		if not (df[column] == np.floor(df[column])).all():
			raise ValueError(f"Non-integer values in {column}.")
		df[column] = df[column].astype(np.int64)

	for column in ("time_sec", "current", "s_rot"):
		df[column] = pd.to_numeric(df[column], errors="raise")
		if not np.isfinite(df[column].to_numpy()).all():
			raise ValueError(f"Non-finite values in {column}.")

	pattern = re.compile(
		r"^(img|ref)_(\d+)_(\d+)\.tiff?$", re.IGNORECASE
	)
	files = {"img": [], "ref": []}
	with os.scandir(folder) as entries:
		for entry in entries:
			match = pattern.fullmatch(entry.name)
			if match and entry.is_file():
				identifier = match.group(1).lower()
				order = (int(match.group(2)), int(match.group(3)))
				files[identifier].append((order, entry.path))

	df["image_file"] = ""
	for identifier, candidates in files.items():
		candidates.sort(key=lambda item: item[0])
		orders = [order for order, _ in candidates]
		if len(set(orders)) != len(orders):
			raise ValueError(
				f"Ambiguous filenames: duplicate numeric {identifier} "
				f"(batch, frame) identifiers in {folder}."
			)

		mask = df["image_identifier"] == identifier
		expected = int(mask.sum())
		if len(candidates) != expected:
			raise ValueError(
				f"{identifier}: {expected} log rows but "
				f"{len(candidates)} matching TIFF files in {folder}. "
				"Cannot safely pair files by acquisition order."
			)
		if expected:
			df.loc[mask, "image_file"] = [
				path for _, path in candidates
			]

	df["image_key"] = df["image_identifier"].map(key_map).astype(np.int64)

	# Follow scanDataset's convention for duplicate/non-increasing times.
	# time_sec retains the original relative timestamps.
	if np.any(np.diff(timestamps) <= 0):
		log.warning(
			"Non-increasing timestamps in %s; normalizing timestamps.",
			logfile,
		)
		timestamps = normalize_timestamps_preserve_increase(timestamps)

	df.index = timestamps
	df["time"] = pd.to_datetime(timestamps, unit="ms")
	if not includeCurrent:
		df.drop(columns=["current"], inplace=True)
	return df
