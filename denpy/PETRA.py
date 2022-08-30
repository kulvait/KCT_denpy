#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 12 18:23:46 2018

Processing of h5 format used to output Petra III data

@author: Vojtech Kulvait

"""
import h5py
import pandas as pd
from bisect import bisect_left


def _insertToDf(df, dat, name):
	time=dat["%s/time"%(name)]
	value=dat["%s/value"%(name)]
	if name == "image_file" and h5py.__version__.startswith("3"):
		value = value.asstr()
	for i in range(len(value)):
		t=time[i]
		v=value[i]
		df.loc[t][name]=v

def beamCurrentDataset(h5file):
	h5 = h5py.File(h5file, 'r')
	timeIndex = filter(lambda t: t!= 0, list(h5["entry/hardware/beam_current/current/time"]))
	df = pd.DataFrame(columns=["time", "current"], index=timeIndex)
	for ind in df.index:
		df.loc[ind]["time"] = pd.to_datetime(ind, unit="ms")
	time = h5["entry/hardware/beam_current/current/time"]
	value = h5["entry/hardware/beam_current/current/value"]
	for i in range(len(value)):
		t = time[i]
		v = value[i]
		if t != 0:
			df.loc[t]["current"] = v 
	return df

def scanDataset(h5file):
	h5 = h5py.File(h5file, 'r')
	data = h5["entry/scan/data"]
	labels = list(data.keys())
	if len(labels) < 1:
		sys.exit("Error: labels count is %d!"%(labels.count))
	df = pd.DataFrame(columns=labels+["current", "time"], index=list(data["%s/time"%(labels[0])][()]))
	for ind in df.index:
		df.loc[ind]["time"] = pd.to_datetime(ind, unit="ms")
	for lab in labels:
		_insertToDf(df, data, lab)
	currentFrame = beamCurrentDataset(h5file)
	for ind in df.index:
		posAfterEq = bisect_left(currentFrame.index, ind)
		if posAfterEq == 0 or posAfterEq == len(df.index) or currentFrame.index[posAfterEq] == ind:
			df.loc[ind]["current"] = currentFrame.iloc[posAfterEq]["current"]
		else:
			t0 = float(currentFrame.index[posAfterEq-1])
			cur0 = float(currentFrame.iloc[posAfterEq-1]["current"])
			t1 = float(currentFrame.index[posAfterEq])
			cur1 = float(currentFrame.iloc[posAfterEq]["current"])
			t = float(ind)
			df.loc[ind]["current"] = cur0 + ((t-t0)/(t1-t0))*(cur1-cur0)
	return df
