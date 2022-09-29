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

#np.arrays are by default in row-major order and they are indexed as follows
#array2d.shape = (dimy, dimx) = (axis0, axis1)
#array3d.shape=(dimz, dimy, dimx)= (axis0, axis1, axis2)


# Get frame from DEN file
# i is the frame index
# Can get a subframe, where row_from is the first row index
# row_to is not included
# col_from is the first col index
# col_to is not included
def getFrame(fileName,
             k,
             row_from=None,
             row_to=None,
             col_from=None,
             col_to=None):
	info = readHeader(fileName)
	columns = info["shape"][-1]
	rows = info["shape"][-2]
	if row_from is None:
		row_from = 0
	if row_to is None:
		row_to = rows
	if col_from is None:
		col_from = 0
	if col_to is None:
		col_to = columns
	f = open(fileName, "rb")
	#For more than 3D arrays
	if info["dimcount"] > 3:
		if len(k) != info["dimcount"] - 2:
			raise TypeError("Index must be a list of dimension %d!" %
			                (info["dimcount"] - 2))
		dim = info["dimspec"][2:]
		blockIncrement = 1
		zindex = 0;
		for i in ind:
			zindex = zindex + k[i]*blockIncrement
			blockIncrement = blockIncrement * dim[i]
	else:
		zindex = k
	offset = info["offset"] + rows * columns * info["elementbytesize"] * zindex
	f.seek(offset, os.SEEK_SET)
	data = np.fromfile(f, dtype=info["type"], count=rows * columns)
	if info["majority"] == 0:
		newdata = data.reshape((rows, columns))
		newdata = newdata[row_from:row_to, col_from:col_to]
	else:
		newdata = data.reshape((columns, rows))
		newdata = np.transpose(newdata)
		newdata = newdata[row_from:row_to, col_from:col_to]
	return (newdata)


def writeFrame(fileName, k, data, force=False):
	if not force and os.path.exists(fileName):
		raise IOError('File %s already exists, no header written' % fileName)
	shape = data.shape
	if len(shape) != 2:
		raise ValueError('Dimension of data should be 2 %d.' % len(shape))
	info = readHeader(fileName)
	columns = info["shape"][-1]
	rows = info["shape"][-2]
	if shape[0] != rows or shape[1] != columns:
		raise ValueError(
		    'There is dimension mismatch between frame (%d, %d) and expected (rows, cols) = (%d, %d) according to header.'
		    % (rows, columns, shape[0], shape[1]))
	#For more than 3D arrays
	if info["dimcount"] > 3:
		if len(k) != info["dimcount"] - 2:
			raise TypeError("Index must be a list of dimension %d!" %
			                (info["dimcount"] - 2))
		k = np.prod(k)
	f = open(fileName, "r+b")
	if info["type"] != data.dtype:
		raise TypeError("Type mismatch between type of DEN %s and type of np.frame %s."%(info["type"], data.dtype))
	data = data.reshape((rows * columns,))
	offset = info["offset"] + rows * columns * info["elementbytesize"] * k
	f.seek(offset, os.SEEK_SET)
	f.write(data.tobytes())
	f.close()


# Indexing of array is [z,y,x] to get x[2] as a third frame
def getNumpyArray(fileName):
	info = readHeader(fileName)
	f = open(fileName, "rb")
	f.seek(info["offset"], os.SEEK_SET)
	if info["majority"] == 0:
		data = np.fromfile(f, dtype=info["type"], count=np.prod(info["shape"]))
		data = data.reshape(info["shape"])
	else:
		data = np.zeros(np.prod(info["shape"]))
		dimspec = info["dimspec"]
		dimx = dimspec[0]
		dimy = dimspec[1]
		dimzextended = np.prod(info["dimspec"][2:])
		data = data.reshape((dimzextended, dimy, dimx))
		for i in range(dimzextended):
			data[i] = getFrame(fileName, i)
	return (data)

def storeNdarrayAsDEN(fileName, dataFrame, ymajor=0, force=False):
	if not force and os.path.exists(fileName):
		raise IOError('File already exists, no data written')
	if not isinstance(dataFrame, np.ndarray):
		raise TypeError('Object dataFrame has to be of type numpy.array')
	dimspec = np.flip(dataFrame.shape, axis=0)
	writeExtendedHeader(
	    fileName,
	    dimspec,
	    elementtype=dataFrame.dtype,
	    ymajor=ymajor,
	    force=force)
	f = open(fileName, "r+b")
	offset = 4096
	f.seek(offset)
	if ymajor == 0:
		f.write(dataFrame.tobytes())
	elif ymajor == 1 and len(dataFrame) > 1:
		if len(dataFrame) == 2:
			f.write(dataFrame.tobytes(order="F"))
		else:
			#Get flat index for better manipulation
			dataFrame.shape = (np.prod(dataFrame.shape[:-2]),
			                   dataFrame.shape[-2], dataFrame.shape[-1])
			frameSize = dataFrame.shape[-2] * dataFrame.shape[-1]
			for k in range(dataFrame.shape[0]):
				newdata = np.array(dataFrame[k])
				seekoffset = offset + frameSize * dataFrame.dtype.itemsize * k
				f.seek(seekoffset, os.SEEK_SET)
				f.write(newdata.tobytes(order="F"))
	else:
		f.close()
		raise TypeError(
		    'Y major order makes no sense for less than two dimensional arrays')
	f.close()
	return True


def writeExtendedHeader(fileName,
                        dimspec,
                        elementtype=np.dtype('<f4'),
                        ymajor=0,
                        force=False):
	if not force and os.path.exists(fileName):
		raise IOError('File %s already exists, no header written' % fileName)
	#Create empty file
	f = open(fileName, "wb")
	f.seek(4095)
	f.write(b"\0")
	header = np.array([
	    0,
	    len(dimspec), elementtype.itemsize, ymajor,
	    npDtypeToDenDataType(elementtype)
	],
	                  dtype='<u2')
	f.seek(0)
	f.write(header.tobytes())
	f.seek(10)
	dimensionSizes = np.array(dimspec, dtype='<u4')
	f.write(dimensionSizes.tobytes())
	f.close()

def writeEmptyDEN(fileName,
                  dimspec,
                  elementtype=np.dtype('<f4'),
                  ymajor=0,
                  force=False):
	if not force and os.path.exists(fileName):
		raise IOError('File %s already exists, no header written' % fileName)
	writeExtendedHeader(fileName, dimspec, elementtype, ymajor, force)
	fileSize = 4096 + np.uint64(np.prod(dimspec)) * elementtype.itemsize
	filesize = np.uint64(fileSize)
	f = open(fileName, "r+b")
	f.seek(int(fileSize - 1))
	f.write(b'\x00')
	f.close()


def denDataTypeToNpDtype(typeID):
	#	 UINT16(2), // 0
	#	 INT16(2), // 1
	#	 UINT32(4), // 2
	#	 INT32(4), // 3
	#	 UINT64(8), // 4
	#	 INT64(8), // 5
	#	 FLOAT32(4), // 6
	#	 FLOAT64(8), // 7
	#	 UINT8(1); // 8
	types = {
	    0: np.dtype('<u2'),
	    1: np.dtype('<i2'),
	    2: np.dtype('<u4'),
	    3: np.dtype('<i4'),
	    4: np.dtype('<u8'),
	    5: np.dtype('<i8'),
	    6: np.dtype('<f4'),
	    7: np.dtype('<f8'),
	    8: np.dtype('<u1')
	}
	return (types[typeID])


def npDtypeToDenDataType(nptype):
	#	 UINT16(2), // 0
	#	 INT16(2), // 1
	#	 UINT32(4), // 2
	#	 INT32(4), // 3
	#	 UINT64(8), // 4
	#	 INT64(8), // 5
	#	 FLOAT32(4), // 6
	#	 FLOAT64(8), // 7
	#	 UINT8(1); // 8
	types = {
	    np.dtype('<u2'): 0,
	    np.dtype('<i2'): 1,
	    np.dtype('<u4'): 2,
	    np.dtype('<i4'): 3,
	    np.dtype('<u8'): 4,
	    np.dtype('<i8'): 5,
	    np.dtype('<f4'): 6,
	    np.dtype('<f8'): 7,
	    np.dtype('<u1'): 8
	}
	return (types[nptype])


def readHeader(fileName):
	header = np.fromfile(fileName, dtype=np.dtype('<u2'), count=3)
	par = {}
	par["h0"] = np.uint32(
	    header[0]).item()  #Convert to default Python integer type
	par["h1"] = np.uint32(header[1]).item(
	)  #see https://stackoverflow.com/questions/9452775/converting-numpy-dtypes-to-native-python-types
	par["h2"] = np.uint32(header[2]).item()
	par["size"] = os.path.getsize(fileName)
	if par["h0"] == 0 and par["size"] > 6:
		header = np.fromfile(fileName, dtype=np.dtype('<u2'), count=2, offset=6)
		par["h3"] = np.uint32(header[0]).item()
		par["h4"] = np.uint32(header[1]).item()
		par["dimcount"] = par["h1"]
		par["elementbytesize"] = par["h2"]
		par["majority"] = par["h3"]
		par["type"] = denDataTypeToNpDtype(par["h4"])
		dimensions = np.fromfile(
		    fileName, dtype=np.dtype('<u4'), count=par["dimcount"], offset=10)
		par["dimspec"] = tuple(dimensions.tolist())
		par["shape"] = tuple(np.flip(dimensions).tolist())
		par["elementscount"] = np.prod(par["shape"])
		par["offset"] = 4096
		par["legacy"] = 0
		if par["elementscount"] * par["elementbytesize"] != par["size"] - par[
		    "offset"]:
			raise TypeError("File %s is not valid DEN file!" % fileName)
	else:
		par["rows"] = np.uint32(
		    header[0]).item()  #Convert to default Python integer type
		par["cols"] = np.uint32(header[1]).item(
		)  #see https://stackoverflow.com/questions/9452775/converting-numpy-dtypes-to-native-python-types
		par["zdim"] = np.uint32(header[2]).item()
		par["shape"] = (par["h2"], par["h0"], par["h1"])  #(dimz, dimy, dimx)
		par["dimspec"] = (par["h1"], par["h0"], par["h2"])  #(dimz, dimy, dimx)
		par["offset"] = 6
		par["dimcount"] = 3
		par["legacy"] = 1
		par["majority"] = 0 #xmajor
		par["elementscount"] = np.prod(par["shape"])
		elementscount = par["elementscount"]
		dataSize = par["size"] - 6
		if dataSize % elementscount != 0:
			raise TypeError("File %s is not valid DEN file!" % fileName)
		if dataSize / elementscount == 2:
			par["type"] = np.dtype('<u2')
			par["elementbytesize"] = 2
		elif dataSize / elementscount == 4:
			par["type"] = np.dtype('<f4')
			par["elementbytesize"] = 4
		elif dataSize / elementscount == 8:
			par["type"] = np.dtype('<f8')
			par["elementbytesize"] = 8
		else:
			raise TypeError("File %s is not valid DEN file!" % fileName)
	if par["elementbytesize"] != par["type"].itemsize:
		raise TypeError(
		    "Byte size of the element %d does not match the byte size of the type %s %d."
		    % (par["elementbytesize"], par["type"].str, par["type"].itemsize))
	elementscount = np.prod(par["shape"])
	return (par)


# Trim frame to the specified dimensions
# Can get a subframe, where row_from is the first row index
# row_to is not included
# col_from is the first col index
# col_to is not included
def trimFrame(frame, row_from, row_to, col_from, col_to):
	newdata = frame[row_from:row_to, col_from:col_to]
	return (newdata)


#Legacy functions
def writeLegacyHeader(fileName, dimx, dimy, dimz, force=False):
	if not force and os.path.exists(fileName):
		raise IOError('File %s already exists, no header written' % fileName)
	dimx = np.uint16(dimx)
	dimy = np.uint16(dimy)
	dimz = np.uint16(dimz)
	outfile = open(fileName, "w")
	header = np.array([dimy, dimx, dimz])
	header = np.array(header, dtype='<u2')
	header.tofile(outfile)
	outfile.close()


def writeEmptyLegacyDEN(fileName, dimx, dimy, dimz, force=False):
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


#Wrong data alignment
def legacyStoreNdarrayAsDEN(fileName, dataFrame, force=False):
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
		raise ValueError('Dimension of dataFrame should be 2 or 3 but is %d.' %
		                 len(dataFrame.shape))
	shape = dataFrame.shape  # Now len is for sure 3
	rows = shape[0]
	columns = shape[1]
	writeDENHeader(
	    fileName, dimx=shape[1], dimy=shape[0], dimz=shape[2], force=force)
	toWrite = dataFrame.astype(np.float32)
	f = open(fileName, "r+b")
	for frame in range(shape[2]):
		newdata = np.array(toWrite[:, :, frame], np.dtype('<f4'))
		newdata = newdata.reshape((rows * columns,))
		# put header in front of image data
		f.seek(6 + rows * columns * frame * 4, os.SEEK_SET)
		f.write(newdata.tobytes())
	f.close()
	return True
