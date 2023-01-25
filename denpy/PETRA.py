#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 12 18:23:46 2018

Processing of h5 format used to output Petra III data

@author: Vojtech Kulvait

"""
import h5py
import pandas as pd
import numpy as np
from bisect import bisect_left


def _insertToDf(df, dat, name):
	time=dat["%s/time"%(name)]
	value=dat["%s/value"%(name)]
	if name == "image_file" and h5py.__version__.startswith("3"):
		value = value.asstr()
	if len(time) != len(value):
		raise IOError("len(time)=%d len(value)=%d for entry %s"%(len(time), len(value), name))
	for i in range(len(value)):
		t=time[i]
		v=value[i]
		if t in df.index:
			df.loc[t][name]=v
		else:
			raise IOError("Can not insert value=%s into the column=%s because related time=%s is not in index"%(v, name, t))

def beamCurrentDataset(h5file):
	h5 = h5py.File(h5file, 'r')
	ID = list(h5["entry/hardware/beam_current/current/time"])
	time = list(pd.to_datetime(i, unit="ms") for i in ID)
	value = list(h5["entry/hardware/beam_current/current/value"])
	agg = list(zip(time, value))
	df = pd.DataFrame(agg, columns=["time", "current"], index = ID)
	#There were duplicates caused problems
	#https://stackoverflow.com/questions/13035764/remove-pandas-rows-with-duplicate-indices
	df = df[~df.index.duplicated(keep='first')]
	#Remove zero index https://stackoverflow.com/questions/13851535/how-to-delete-rows-from-a-pandas-dataframe-based-on-a-conditional-expression
	df.drop(df[df.index == 0].index, inplace=True)
	#Need to sort it for interpolation to work
	df.sort_index(inplace=True)
	return df

def scanDataset(h5file):
	h5 = h5py.File(h5file, 'r')
	data = h5["entry/scan/data"]
	labels = list(data.keys())
	if len(labels) < 1:
		sys.exit("Error: labels count is %d!"%(labels.count))
	#There always shall be image_key entry
	df = pd.DataFrame(columns=labels+["current", "time"], index=list(data["image_key/time"]))
	for ind in df.index:
		df.loc[ind, "time"] = pd.to_datetime(ind, unit="ms")
	for lab in labels:
		try:
			_insertToDf(df, data, lab)
		except IOError:
			print("IOError processing %s"%h5file)
			raise
	currentFrame = beamCurrentDataset(h5file)
	for ind in df.index:
		posAfterEq = bisect_left(currentFrame.index, ind)
		if posAfterEq == len(currentFrame.index):
			df.loc[ind]["current"] = currentFrame.iloc[posAfterEq-1]["current"]
		elif posAfterEq == 0 or currentFrame.index[posAfterEq] == ind:
			df.loc[ind]["current"] = currentFrame.iloc[posAfterEq]["current"]
		else:
			t0 = float(currentFrame.index[posAfterEq-1])
			cur0 = float(currentFrame.iloc[posAfterEq-1]["current"])
			t1 = float(currentFrame.index[posAfterEq])
			cur1 = float(currentFrame.iloc[posAfterEq]["current"])
			t = float(ind)
			cur = cur0 + ((t-t0)/(t1-t0))*(cur1-cur0)
			if np.isnan(cur):
				print("Cur is NaN t0=%f cur0=%f t1=%f cur1=%f t=%f ind=%d posAfterEq=%d len(currentFrame)=%d!"%(t0, cur0, t1, cur1, t, ind, posAfterEq, len(currentFrame)))
				raise ValueError("Interpolation error producing NaN beam curent.")
			df.loc[ind]["current"] = cur
	return df
