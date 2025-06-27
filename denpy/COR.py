#!/usr/bin/env python
#-*- coding: utf-8 -*-
"""
Created 2024-2025

Implementation of various techniques for center of rotation detection in parallel ray geometry

@author: Vojtech Kulvait
@license: GNU General Public License v3.0 (GPLv3)
"""
import os
import numpy as np
import logging
import sys
import scipy
from scipy.ndimage import gaussian_filter
from scipy.signal import savgol_filter
import matplotlib.pyplot as plt
import matplotlib.ticker as plticker


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

#Center of rotation is estimated as a offset from center 0.5xdim so that COR = 0.5xdim + offset
#For convenience x dimension is measured from the edge of the 0th pixel so that the
#x coordinate of the pixel center is x(i) = i + 0.5
#It works with sinograms, where(ydim, xdim) = shape.sinogram and
#ydim... measurements for different angles accessible by sinogram[j]
#xdim... measurements for individual detector pixels

def maskLongerTail(sin, offset):
	"""
	Masks the tails of a sinogram based on a given offset from the center.

	When the center of rotation estimate is offset from the middle of the sinogram,
	the tails on one side become "unbalanced" (longer). This function masks (zeros)
	the longer tail part to balance the sinogram data for improved center estimation.

	Assumptions:
	- The offset is expected to be rounded to the nearest multiple of 0.5 pixels.
	- The sinogram data `sin` is a 2D ndarray with shape (angles, xdim) or 1D ndarray of length xdim.
	- Pixel indices start at 0, with pixel centers spaced by 1 unit.

	Parameters:
	----------
	sin : np.ndarray
		The input sinogram data array.
		Shape can be (ydim, xdim) or (xdim,) for single angle.

	offset : float
		The estimated shift from the sinogram center (in pixels).
		Positive values mean the center is shifted right, negative left.

	Returns:
	-------
	np.ndarray
		A copy of the input sinogram with the longer tail masked (set to zero).
		The shape is preserved.

	Notes:
	-----
	- The length of the masked tail is `2 * abs(offset)` pixels,
	  reflecting that a 0.5 pixel offset corresponds to masking 1 pixel.
	- If offset is zero, the input sinogram is returned unchanged.
	"""
	if not isinstance(sin, np.ndarray):
		raise ValueError(
			"Object sin shall be ndarray but it is %s" % (type(sin)))
	# Round offset to nearest 0.5 multiple
	offset = round(offset * 2) / 2
	# Do nothing on zero offset
	if offset == 0:
		return sin
	sin = sin.copy()
	reshape = False
	if len(sin.shape) == 1:
		xdim = sin.shape[0]
		sin.shape = (1, xdim)
		reshape = True
	else:
		xdim = sin.shape[1]
	# When offset is multiple of 0.5 then masked length is 2*abs(offset)
	maskedLen = int(2 * np.abs(offset))
	if offset < 0:
		sin[:, -maskedLen:] = 0
	else:
		sin[:, :maskedLen] = 0
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
	#Truncation towards 0 is safer on the boundary
	initShift = int(2 * initOffset)
	if initShift >= xdim or initShift <= -xdim:
		raise ValueError(f"For input parameters xdim={xdim}, initOffset={initOffset}, searchDiameter={searchDiameter}, balanced={balanced}:\n"
			f"invalid value of initShift=int(2*initOffset)={initShift} is out of range (-xdim, xdim).\n" )
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
		log.info("Iteration %d with xdim=%d, init_offset=%0.1f, nrmord=%d, search_diameter=%s, balanced=%s, peakSharpness=%0.2f" %
				 (iteration, xdim, init_offset, nrmord, search_diameter, balanced, peakSharpness))
		maskingParams = computeMaskingParameters(xdim,
												 init_offset,
												 searchDiameter=search_diameter,
												 balanced=balanced)
		log.info("Searching diameter %d from init_offset=%0.1f in xdim=%d maskCount=[%d, %d]" %
					 (maskingParams["searchDiameter"], maskingParams["initOffset"], xdim, maskingParams["maskLeft"],
					  maskingParams["maskRight"]))
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
	# --- Warning if minimizer is at the boundary of the search range
	if shiftEstimate_f == xdim - 1 or shiftEstimate_f == 0:
		# This is a warning, not an error, so we use print instead of raise
		fname = os.path.basename(__file__) if '__file__' in globals() else "unknown"
		print(
			f"WARNING denpy.COR.sinogram_consistency_detection360 ({fname}): "
			f"Minimizer lies on boundary of global search range "
			f"(shiftEstimate_f={shiftEstimate_f}, full xdim={xdim}, local IND=({IND[0]}, {IND[-1]}))",
			file=sys.stderr
		)
	if verbose:
		#plt.title(
		#	 "Init offset = %0.1f in green first half integer estimate %0.1f in red."
		#	 % (init_offset, 0.5 * shiftEstimate_cur))
		plt.title(
			"Initial unbinned offset estimate %0.1f in red."
			% (0.5 * shiftEstimate_cur))
		plt.plot(0.5*IND, normedValues)
		#plt.axvline(x=2 * init_offset, color="green", linewidth=1)
		plt.axvline(x=shiftEstimate_cur*0.5, color="red", linewidth=1)
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
		#plt.title(
		#	 "Offset estimate first %0.1f in green and converged %0.1f in red." %
		#	 (init_offset, 0.5 * shiftEstimate_cur))
		plt.gca().xaxis.set_major_locator(plt.MultipleLocator(2))
		plt.gca().xaxis.set_minor_locator(plt.MultipleLocator(1))
		plt.title(
			"Unbinned offset estimate %0.1f in red." %
			(0.5 * shiftEstimate_cur))
		plt.plot(0.5*IND, normedValues)
		#plt.axvline(x=shiftEstimate_f, color="green", linewidth=1)
		plt.axvline(x=shiftEstimate_cur*0.5, color="red", linewidth=1)
		plt.xlim(shiftEstimate_cur*0.5 - 10, shiftEstimate_cur*0.5 + 10)
		plt.ylabel("Minimizer value")
		plt.xlabel("Shift offset")
		plt.show()
	#Try 1D interpolation with splines
	#normedValues are values of quantity to minimize on integer points IND
	interpolation = scipy.interpolate.interp1d(IND, normedValues, kind='cubic')
	#Now subpixel precision
	try:
		searchDiameter = 5
		maskingParams = computeMaskingParameters(xdim, shiftEstimate_cur * 0.5,
												 searchDiameter, balanced)
		IND_ORIG = IND
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
		plt.gca().xaxis.set_major_locator(plt.MultipleLocator(2))
		plt.gca().xaxis.set_minor_locator(plt.MultipleLocator(1))
		plt.title(
			"Unbinned offset estimate interpolation %0.2f in red." %
			(ESTIMATE_INTERP))
		plt.plot(0.5*IND_ORIG, normedValues, label="Original")
		plt.plot(0.5*IND, interpValues, label="Interpolation")
		plt.axvline(x=ESTIMATE_INTERP, color="red", linewidth=1)
		plt.xlim(shiftEstimate_cur*0.5 - 10, shiftEstimate_cur*0.5 + 10)
		plt.legend()
		plt.show()
	return (ESTIMATE_OFFSET, ESTIMATE_INTERP, convergingSequence,
			minimizerValue, peakSharpnessGlobal, peakSharpnessLocal, IND, interpValues)
