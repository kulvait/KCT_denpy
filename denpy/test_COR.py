#!/usr/bin/env python
#-*- coding: utf-8 -*-"
"""
Created 2024

Test of Implementation of various techniques for center of rotation detection in parallel ray geometry

@author: Vojtech Kulvait
"""
import unittest
import numpy as np

from denpy import COR
from denpy.COR import maskLongerTail
from denpy.COR import computeMaskingParameters


def checkMaskingParameters(par):
	xdim = par["xdim"]
	assert (xdim > 0)
	initOffset = par["initOffset"]
	searchDiameter = par["searchDiameter"]
	assert (searchDiameter > 0)
	assert (np.abs(initOffset * 2) + searchDiameter < xdim)
	sliceMask = par["slice"]
	balanced = par["balanced"]
	A = COR.getShiftArray(xdim, initOffset)
	#	print("A=%s"%(A))
	if balanced:  #When searching distances shall be balanced
		A_masked = A[sliceMask]
		#		print("A_%s = %s"%(sliceMask, A_masked))
		assert (-A_masked.min() == A_masked.max())
	IND = par["ind"]
	maskLeft = par["maskLeft"]
	maskRight = par["maskRight"]
	assert (maskRight + maskLeft < xdim)
	#Now create mask by shifting and test wether it corresponds to created mask
	initShift = int(round(initOffset * 2, 0))
	AF = np.flip(A)
	AFL = np.roll(AF, initShift + searchDiameter)
	#	print("searchDiameter=%s A+AFL=%s"%(searchDiameter, A+AFL))
	ALPHALEFT = ((A + AFL) == searchDiameter)
	AFR = np.roll(AF, initShift - searchDiameter)
	ALPHARIGHT = ((A + AFR) == -searchDiameter)
	ALPHA = np.multiply(ALPHALEFT, ALPHARIGHT)
	#	print("ALPHALEFT=%s"%(ALPHALEFT))
	#	print("ALPHARIGHT=%s"%(ALPHARIGHT))
	#	print("ALPHA=%s"%(ALPHA))
	#ALPHA is now maximal mask False = masked, True = not masked
	alphaLeft = ALPHA.tolist().index(1)
	alphaRight = np.flip(ALPHA).tolist().index(1)
	#	print("alphaLeft=%d maskLeft=%d"%(alphaLeft, maskLeft))
	#	print("alphaLeft=%d maskLeft=%d"%(alphaRight, maskRight))
	if not balanced:
		assert (alphaLeft == maskLeft and alphaRight == maskRight)
	else:
		assert (alphaLeft == maskLeft or alphaRight == maskRight)
		assert (alphaLeft <= maskLeft and alphaRight <= maskRight)


class CORTest(unittest.TestCase):

	def test_maskLongerTail(self):
		arrayA = np.array([1, 2, 3, 4, 5])
		arrayB = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
		arrayB.shape = (1, 10)
		np.testing.assert_array_equal(maskLongerTail(arrayA, -0.5),
		                              np.array([1, 2, 3, 4, 0]))
		np.testing.assert_array_equal(maskLongerTail(arrayA, 0.5),
		                              np.array([0, 2, 3, 4, 5]))
		np.testing.assert_array_equal(maskLongerTail(arrayA, -1.2),
		                              np.array([1, 2, 3, 0, 0]))
		np.testing.assert_array_equal(maskLongerTail(arrayA, 1),
		                              np.array([0, 0, 3, 4, 5]))
		np.testing.assert_array_equal(maskLongerTail(arrayA, 2.5),
		                              np.array([0, 0, 0, 0, 0]))
		np.testing.assert_array_equal(
		    maskLongerTail(arrayB, -0.5),
		    np.array([[1, 2, 3, 4, 5, 6, 7, 8, 9, 0]]))
		np.testing.assert_array_equal(
		    maskLongerTail(arrayB, 0.5),
		    np.array([[0, 2, 3, 4, 5, 6, 7, 8, 9, 10]]))
		np.testing.assert_array_equal(
		    maskLongerTail(arrayB, -1.2),
		    np.array([[1, 2, 3, 4, 5, 6, 7, 8, 0, 0]]))
		np.testing.assert_array_equal(
		    maskLongerTail(arrayB, 1),
		    np.array([[0, 0, 3, 4, 5, 6, 7, 8, 9, 10]]))
		np.testing.assert_array_equal(
		    maskLongerTail(arrayB, 2.5),
		    np.array([[0, 0, 0, 0, 0, 6, 7, 8, 9, 10]]))
		np.testing.assert_array_equal(
		    maskLongerTail(arrayB, 0.01),
		    np.array([[1, 2, 3, 4, 5, 6, 7, 8, 9, 10]]))

	def test_computeMaskingParameters(self):
		checkMaskingParameters(computeMaskingParameters(100, -5.0, 10, False))
		checkMaskingParameters(computeMaskingParameters(100, -5.0, 10, True))
		checkMaskingParameters(computeMaskingParameters(99, -5.0, 10, False))
		checkMaskingParameters(computeMaskingParameters(99, 5.5, 10, True))
		checkMaskingParameters(
		    computeMaskingParameters(1099, -5.0, balanced=False))
		checkMaskingParameters(
		    computeMaskingParameters(1100, -5.0, balanced=True))


if __name__ == '__main__':
	unittest.main()
