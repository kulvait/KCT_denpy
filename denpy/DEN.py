#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 12 18:23:46 2018

DEN and DICOM IO manipulation

@author: Vojtech Kulvait
@author: Enno Schieferecke

"""
import numpy as np
import os


# Get frame from DEN file
# i is the frame index
# Can get a subframe, where row_from is the first row index
# row_to is not included
# col_from is the first col index
# col_to is not included
def getFrame(fileName, i, row_from=None, row_to=None, col_from=None, col_to=None):
	header = np.fromfile(fileName, np.dtype('<i2'), 3)
	rows = np.uint32(header[0])
	columns = np.uint32(header[1])
	if row_from is None:
		row_from = 0
	if row_to is None:
		row_to = rows
	if col_from is None:
		col_from = 0
	if col_to is None:
		col_to = columns
	f = open(fileName, "rb")
	f.seek(6 + rows * columns * 4 * i, os.SEEK_SET)
	data = np.fromfile(f, np.dtype('<f4'), rows * columns)
	newdata = data.reshape((rows, columns))
	newdata = newdata[row_from:row_to, col_from:col_to]
	return(newdata)

# Indexing of array is [z,y,x] to get x[2] as a third frame


def getNumpyArray(fileName):
	header = readHeader(fileName)
	f = open(fileName, "rb")
	f.seek(6, os.SEEK_SET)
	data = np.fromfile(f, np.dtype('<f4'), header[
					   "rows"] * header["cols"] * header["zdim"])
	data = data.reshape((header["zdim"], header["rows"], header["cols"]))
	return(data)


def storeNdarrayAsDEN(fileName, dataFrame, force=False):
	if not force and os.path.exists(fileName):
		raise IOError('File already exists, no data written')
	if not isinstance(dataFrame, np.ndarray):
		raise TypeError('Object dataFrame has to be of type numpy.array')
	if len(dataFrame.shape) == 1:
		print('Dimension = 1, expected >= 2')
		return False
	elif len(dataFrame.shape) == 2:
		dataFrame = np.expand_dims(dataFrame, axis=2)
	elif len(dataFrame.shape) > 3:
		raise ValueError(
			'Dimension of dataFrame should be 2 or 3 but is %d.' % len(dataFrame.shape))
	shape = dataFrame.shape  # Now len is for sure 3
	rows = shape[0]
	columns = shape[1]
	writeDENHeader(fileName, dimx=shape[1], dimy=shape[0], dimz=shape[2],
				  force=force)
	toWrite = dataFrame.astype(np.float32)
	f = open(fileName, "r+b")
	for frame in range(shape[2]):
		newdata = np.array(toWrite[:, :, frame], np.dtype('<f4'))
		newdata = newdata.reshape((rows * columns, ))
# put header in front of image data
		f.seek(6 + rows * columns * frame * 4, os.SEEK_SET)
		f.write(newdata.tobytes())
	f.close()
	return True


def storeNdarrayAsFloatDEN(fileName, dataFrame, force=False):
	if not force and os.path.exists(fileName):
		raise IOError('File already exists, no data written')
	if not isinstance(dataFrame, np.ndarray):
		raise TypeError('Object dataFrame has to be of type numpy.array')
	if len(dataFrame.shape) == 1:
		print('Dimension = 1, expected >= 2')
		return False
	elif len(dataFrame.shape) == 2:
		dataFrame = np.expand_dims(dataFrame, axis=2)
	elif len(dataFrame.shape) > 3:
		raise ValueError(
			'Dimension of dataFrame should be 2 or 3 but is %d.' % len(dataFrame.shape))
	shape = dataFrame.shape  # Now len is for sure 3
	rows = shape[0]
	columns = shape[1]
	writeDENHeader(fileName, dimx=shape[1], dimy=shape[0], dimz=shape[2],
				  force=force)
	toWrite = dataFrame.astype(np.float32)
	f = open(fileName, "r+b")
	for frame in range(shape[2]):
		newdata = np.array(toWrite[:, :, frame], np.dtype('<f4'))
		newdata = newdata.reshape((rows * columns, ))
# put header in front of image data
		f.seek(6 + rows * columns * frame * 4, os.SEEK_SET)
		f.write(newdata.tobytes())
	f.close()
	return True

def storeNdarrayAsDoubleDEN(fileName, dataFrame, force=False):
	if not force and os.path.exists(fileName):
		raise IOError('File already exists, no data written')
	if not isinstance(dataFrame, np.ndarray):
		raise TypeError('Object dataFrame has to be of type numpy.array')
	if len(dataFrame.shape) == 1:
		print('Dimension = 1, expected >= 2')
		return False
	elif len(dataFrame.shape) == 2:
		dataFrame = np.expand_dims(dataFrame, axis=2)
	elif len(dataFrame.shape) > 3:
		raise ValueError(
			'Dimension of dataFrame should be 2 or 3 but is %d.' % len(dataFrame.shape))
	shape = dataFrame.shape  # Now len is for sure 3
	rows = shape[0]
	columns = shape[1]
	writeDENHeader(fileName, dimx=shape[1], dimy=shape[0], dimz=shape[2],
				  force=force)
	toWrite = dataFrame.astype(np.float64)
	f = open(fileName, "r+b")
	for frame in range(shape[2]):
		newdata = np.array(toWrite[:, :, frame], np.dtype('<f8'))
		newdata = newdata.reshape((rows * columns, ))
# put header in front of image data
		f.seek(6 + rows * columns * frame * 8, os.SEEK_SET)
		f.write(newdata.tobytes())
	f.close()
	return True


def writeFrame(fileName, k, data, force=False):
	if not force and os.path.exists(fileName):
		raise IOError('File %s already exists, no header written' % fileName)
	shape = data.shape
	if len(shape) != 2:
		raise ValueError('Dimension of data should be 2 %d.' % len(shape))
	header = np.fromfile(fileName, np.dtype('<i2'), 3)
	rows = np.uint64(header[0])
	columns = np.uint64(header[1])
	k = np.uint64(k)
	if shape[0] != rows or shape[1] != columns:
		raise ValueError('There is dimension mismatch between frame (%d, %d) and expected (rows, cols) = (%d, %d) according to header.' %
						 (rows, columns, shape[0], shape[1]))
	f = open(fileName, "r+b")
	data = np.array(data, np.dtype('<f4'))
	data = data.reshape((rows * columns, ))
	f.seek(6 + rows * columns * k * 4, os.SEEK_SET)
	data.tofile(f)
	f.close()


def writeDENHeader(fileName, dimx, dimy, dimz, force=False):
	if not force and os.path.exists(fileName):
		raise IOError('File %s already exists, no header written' % fileName)
		dimx = np.uint(dimx)
		dimy = np.uint(dimy)
		dimz = np.uint(dimz)
	outfile = open(fileName, "w")
	header = np.array([dimy, dimx, dimz])
	header = np.array(header, dtype='<i2')
	header.tofile(outfile)
	outfile.close()

def writeEmptyDEN(fileName, dimx, dimy, dimz, force=False):
	if not force and os.path.exists(fileName):
		raise IOError('File %s already exists, no header written' % fileName)
		dimx = np.uint(dimx)
		dimy = np.uint(dimy)
		dimz = np.uint(dimz)
	outfile = open(fileName, "w")
	header = np.array([dimy, dimx, dimz])
	header = np.array(header, dtype='<i2')
	header.tofile(outfile)
	fileSize = dimx * dimy * dimz * 4 + 6
	outfile.seek(fileSize - 1)
	outfile.write('\0')
	outfile.close()


def readHeader(file):
	header = np.fromfile(file, np.dtype('<i2'), 3)
	par = {}
	par["rows"] = np.uint32(header[0])
	par["cols"] = np.uint32(header[1])
	par["zdim"] = np.uint32(header[2])
	return(par)

# Trim frame to the specified dimensions
# Can get a subframe, where row_from is the first row index
# row_to is not included
# col_from is the first col index
# col_to is not included


def trimFrame(frame, row_from, row_to, col_from, col_to):
	newdata = frame[row_from:row_to, col_from:col_to]
	return(newdata)
