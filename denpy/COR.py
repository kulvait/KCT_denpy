#!/usr/bin/env python
#-*- coding: utf-8 -*-
"""
Created 2024

Implementation of various techniques for center of rotation detection in parallel ray geometry

@author: Vojtech Kulvait
"""
import numpy as np
import logging
import sys
import scipy
from scipy.ndimage import gaussian_filter
from scipy.signal import savgol_filter
import matplotlib.pyplot as plt


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

#Center of rotation is estimated as a offset from center 0.5xdim so that COR = 0.5xdim + offset
#For convenience x dimension is measured from the edge of the 0th pixel so that the
#x coordinate of the pixel center is x(i) = i + 0.5
#It works with sinograms, where(ydim, xdim) = shape.sinogram and
#ydim... measurements for different angles accessible by sinogram[j]
#xdim... measurements for individual detector pixels


#When search start is offsetted from the center the tails shall be balanced
#Make zeros at the positions where there is no matching tail
#Rounds shifts to multiplies of 0.5
#Lets expect offset is an integer multiple of 0.5, then 0.5*xdim + offset is the position in the array that will be ballanced
def maskLongerTail(sin, offset):
	if not isinstance(sin, np.ndarray):
		raise ValueError(
		    "Object sin shall be ndarray but it is %s" % (type(sin)))
	sin = sin.copy()
	reshape = False
	if len(sin.shape) == 1:
		xdim = sin.shape[0]
		sin.shape = (1, xdim)
		reshape = True
	else:
		xdim = sin.shape[1]
	midpoint = 0.5 * xdim
	centerPoint = midpoint + offset
	center_two = round(centerPoint * 2, 0)
	betweenPixels = True  #If midpoint lies just between the pixels
	if center_two % 2 == 1:
		betweenPixels = False
	centerPoint = center_two / 2
	centerPixel = int(center_two // 2)
	if centerPoint == midpoint:
		return sin
	elif centerPoint < midpoint:  #Right tail will be censored
		if betweenPixels:  #centerpoint is centerPixel+0.5
			censorIndex = centerPixel * 2
		else:
			censorIndex = 1 + centerPixel * 2
		sin[:, censorIndex:] = 0
	else:  #Left tail will be censored
		if betweenPixels:  #centerpoint is centerPixel+0.5
			censorLength = 2 * centerPixel - xdim
		else:
			censorLength = 2 * centerPixel + 1 - xdim
		sin[:, :censorLength] = 0
	if reshape:
		sin.shape = (xdim)
	return sin


#Vizualize offset based coordinate in array in the middle of each pixel
def getShiftArray(xdim, offset):
	A = np.arange(xdim)
	leftValue = 0.5 * (1 - xdim)
	A = np.arange(xdim)
	A = A + leftValue - offset
	return A


def sliceFromOffsets(xdim, maskLeft=0, maskRight=0):
	if type(maskLeft) != int or type(maskRight) != int or type(xdim) != int:
		ValueError(
		    "Type of xdim=%s maskleft=%s and maskRight=%s shall be integer!" %
		    (xdim, maskLeft, maskRight))
	sliceLeft = None
	sliceRight = None
	if maskLeft != 0:
		sliceLeft = int(maskLeft)
	if maskRight != 0:
		sliceRight = xdim - maskRight
	return slice(sliceLeft, sliceRight)


#masking parameters for check alignment methods to be ballanced checking same number of pixels in each step
#xdim is the length of array to search
#initOffset shall be multiple of 0.5 otherwise it will be rounded
#searchDiameter number of pixels around initOffset to search shall be multiple of 1, otherwise will be rounded
#balanced whether shall the search area be ballanced around initOffset, causes greeder masking
def computeMaskingParameters(xdim,
                             initOffset=0.0,
                             searchDiameter=None,
                             balanced=False):
	initShift = int(round(initOffset * 2, 0))
	if initShift >= xdim or initShift <= -xdim:
		raise ValueError(
		    "Impossible 2*initOffset=initShift<xdim %d>=%d" % (initShift, xdim))
	#How many values to mask on left and right
	maskLeft = max(0, initShift)
	maskRight = max(0, -initShift)
	#What is maximum and minimum shift to still have at least one cell to compare
	maxShift = xdim - 1
	minShift = -maxShift
	#with respect to initShift what are those
	maxShift_relative = maxShift - initShift
	minShift_relative = minShift - initShift
	#These would incude large shifts and that is undesirable
	maxDiameter = min(maxShift_relative, -minShift_relative)
	if type(searchDiameter) == float:
		searchDiameter = int(searchDiameter * maxDiameter)
	if searchDiameter is None:
		searchDiameter = int(0.33 * maxDiameter)
		#Arbitrary choice to have some search area left while still large search area
		if searchDiameter < 1:
			raise ValueError(
			    "Estimated search diameter is too small for given values!")
	#print("initShift=%d maxShift_relative=%d minShift_relative=%d maxDiameter=%d searchDiameter=%d"%(initShift, maxShift_relative, minShift_relative, maxDiameter, searchDiameter))
	maskLeft = int(max(0, initShift + searchDiameter))
	maskRight = int(max(0, -initShift + searchDiameter))
	initPos = 0.5 * xdim + initOffset
	if balanced:
		#We need to initPos-maskLeft == xdim -maskRight - initPos
		if initPos - maskLeft < xdim - maskRight - initPos:
			maskRight = int(xdim - 2 * initPos + maskLeft)
		else:
			maskLeft = int(2 * initPos + maskRight - xdim)
	if maskLeft > initPos or maskRight > (
	    xdim - initPos) or maskLeft + maskRight >= xdim:
		raise ValueError("Impossible masking conditions")
	IND = np.arange(initShift - searchDiameter, initShift + searchDiameter + 1)
	sliceMask = sliceFromOffsets(xdim, maskLeft, maskRight)
	#Can be used to  slice https://docs.python.org/3/library/functions.html#slice
	par = {}
	par["xdim"] = xdim
	par["initOffset"] = initOffset
	par["searchDiameter"] = searchDiameter
	par["balanced"] = balanced
	par["ind"] = IND
	par["maskLeft"] = maskLeft
	par["maskRight"] = maskRight
	par["slice"] = sliceMask
	return par


def monotonicityChanges(convergingSequence):
	#First check that the converging sequence does not have max and min ARGMIN
	vals = [x[1] for x in convergingSequence]
	if len(vals) < 3:
		return 0
	nondecreasing = True
	if vals[1] < vals[0]:
		nondecreasing = False
	monotonicityChanges = 0
	oldval = vals[1]
	for v in vals[2:]:
		if nondecreasing:
			if v < oldval:
				nondecreasing = False
				monotonicityChanges += 1
		else:
			if v >= oldval:
				noncecreasing = True
				monotonicityChanges += 1
		oldval = v
	return monotonicityChanges


def running_mean(x, N):
	out = np.zeros_like(x, dtype=np.float64)
	dim_len = x.shape[0]
	for i in range(dim_len):
		if N % 2 == 0:
			a, b = i - (N - 1) // 2, i + (N - 1) // 2 + 2
		else:
			a, b = i - (N - 1) // 2, i + (N - 1) // 2 + 1
		#cap indices to min and max indices
		a = max(0, a)
		b = min(dim_len, b)
		out[i] = np.mean(x[a:b])
	return out


def searchShiftRange(sa, sb, IND, maskSlice, nrmord=0):
	normedValues = np.zeros(len(IND))
	x1 = sa[:, maskSlice]
	x1 = x1.flatten()
	for i in range(len(IND)):
		x2 = np.roll(sb, IND[i], axis=1)[:, maskSlice]
		x2 = x2.flatten()  #For corrcoef computation
		if nrmord == 0:
			normedValues[i] = np.corrcoef(x1, x2)[0, 1]
		else:
			normedValues[i] = np.linalg.norm(x1 + x2, ord=nrmord)
	ARGMIN = np.argmin(normedValues)
	shiftEstimate = IND[ARGMIN]
	return shiftEstimate, ARGMIN, normedValues


def measurePeakSharpness(normedValues, ARGMIN):
	normedValues = normedValues - normedValues.mean()
	normedValuesMean20 = running_mean(normedValues, 20)
	normedValues20 = normedValues - normedValuesMean20
	std = np.std(normedValues)
	peakSharpnessGlobal = -normedValues[ARGMIN] / std
	peakSharpnessLocal = -normedValues20[ARGMIN] / std
	return peakSharpnessGlobal, peakSharpnessLocal


def sinogram_consistency_detection360(sinogram,
                                      init_offset=0.0,
                                      nrmord=0,
                                      search_diameter=None,
                                      balanced=False,
                                      verbose=False):
	xdim = sinogram.shape[1]
	halfAngleCount = sinogram.shape[0] // 2
	sa = sinogram[0:halfAngleCount]
	sb = -np.flip(sinogram[halfAngleCount:], axis=1)
	peakSharpness = 0
	iteration = 0
	while iteration < 6 and peakSharpness < 5.0:
		maskingParams = computeMaskingParameters(xdim,
		                                         init_offset,
		                                         searchDiameter=search_diameter,
		                                         balanced=balanced)
		#	log.info("Searching diameter %d from init_offset=%0.1f in xdim=%d maskCount=[%d, %d]" %
		#	         (maskingParams["searchDiameter"], maskingParams["initOffset"], xdim, maskingParams["maskLeft"],
		#	          maskingParams["maskRight"]))
		IND = maskingParams["ind"]
		maskSlice = maskingParams["slice"]
		shiftEstimate_old = 0
		shiftEstimate_cur, ARGMIN, normedValues = searchShiftRange(
		    sa, sb, IND, maskSlice, nrmord)
		peakSharpnessGlobal, peakSharpnessLocal = measurePeakSharpness(
		    normedValues, ARGMIN)
		#print("shiftEstimate_cur=%f peakSharpnessLocal=%f peakSharpnessGlobal=%f"%(shiftEstimate_cur, peakSharpnessLocal, peakSharpnessGlobal))
		peakSharpness = peakSharpnessGlobal
		init_offset = 0.5 * shiftEstimate_cur
		iteration = iteration + 1
	shiftEstimate_f = IND[ARGMIN]
	if verbose:
		plt.title(
		    "Init offset = %0.1f in green first half integer estimate %0.1f in red."
		    % (init_offset, 0.5 * shiftEstimate_cur))
		plt.plot(IND, normedValues)
		plt.axvline(x=2 * init_offset, color="green", linewidth=1)
		plt.axvline(x=shiftEstimate_cur, color="red", linewidth=1)
		plt.ylabel("Minimizer value")
		plt.xlabel("Shift offset")
		plt.show()
	#Find noninterpolated minimum
	searchDiameter = 100
	IND = maskingParams["ind"]
	convergingSequence = []  #Index to array IND and value
	minimizerValue = 1e10
	shiftEstimate_old = -shiftEstimate_cur
	while len(convergingSequence) < 1 or (
	    shiftEstimate_cur != shiftEstimate_old and
	    monotonicityChanges(convergingSequence) < 3):
		try:
			shiftEstimate_old = shiftEstimate_cur
			maskingParams = computeMaskingParameters(
			    xdim, shiftEstimate_cur * 0.5, searchDiameter, balanced)
			IND = maskingParams["ind"]
			maskSlice = maskingParams["slice"]
			shiftEstimate_cur, ARGMIN, normedValues = searchShiftRange(
			    sa, sb, IND, maskSlice, nrmord)
			peakSharpnessGlobal, peakSharpnessLocal = measurePeakSharpness(
			    normedValues, ARGMIN)
			convergingSequence.append(
			    (ARGMIN, shiftEstimate_cur,
			     peakSharpnessGlobal + peakSharpnessLocal))
			#print("shiftEstimate_cur=%f peakSharpnessLocal=%f peakSharpnessGlobal=%f"%(shiftEstimate_cur, peakSharpnessLocal, peakSharpnessGlobal))
		except:
			ESTIMATE_OFFSET = None
			break
	else:  #On no break
		ESTIMATE_OFFSET = 0.5 * shiftEstimate_cur
	peakSharpnessGlobal, peakSharpnessLocal = measurePeakSharpness(
	    normedValues, ARGMIN)
	if verbose:
		plt.title(
		    "Offset estimate first %0.1f in green and converged %0.1f in red." %
		    (init_offset, 0.5 * shiftEstimate_cur))
		plt.plot(IND, normedValues)
		plt.axvline(x=shiftEstimate_f, color="green", linewidth=1)
		plt.axvline(x=shiftEstimate_cur, color="red", linewidth=1)
		plt.ylabel("Minimizer value")
		plt.xlabel("Shift offset")
		plt.show()
	#Try 1D interpolation with splines
	interpolation = scipy.interpolate.interp1d(IND, normedValues, kind='cubic')
	#Now subpixel precision
	try:
		searchDiameter = 5
		maskingParams = computeMaskingParameters(xdim, shiftEstimate_cur * 0.5,
		                                         searchDiameter, balanced)
		IND = maskingParams["ind"]
		IND = np.linspace(IND[0], IND[-1], 1001)
		#We intentionally keep previous maskSlice to avoid interpolation errors so not	maskSlice = par["slice"]
		splineOrder = 3
		#Procedure to compute normedValues for given sequence of IND
		interpValues = np.zeros(len(IND))
		for i in range(len(IND)):
			interpValues[i] = interpolation(IND[i])
		ARGMIN = np.argmin(interpValues)
		ESTIMATE_INTERP = 0.5 * IND[ARGMIN]
	except Exception as ex:
		ESTIMATE_INTERP = None
		raise
	#plt.axvline(x=0.5*MAX -  0.5* IND[ARGMII], color="blue", linewidth=3)
	#plt.show()
	if verbose:
		plt.title("ESTIMATE_SUBPIX and ESTIMATE_INTERP")
		plt.plot(IND, interpValues, label="Interpolation")
		plt.legend()
		plt.show()
	return (ESTIMATE_OFFSET, ESTIMATE_INTERP, convergingSequence,
	        minimizerValue, peakSharpnessGlobal, peakSharpnessLocal)
