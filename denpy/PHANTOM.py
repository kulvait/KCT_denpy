#!/usr/bin/env python
#-*- coding: utf-8 -*-
"""
Creating Shepp-Logan like phantoms in 2D and 3D storing them as DEN
@author: Vojtěch Kulvait
(c) 2023

Original inspiration for the code is
https://gist.github.com/blink1073/6e417c726b0dd03d5ea0
Code was aparently inspirred by the contribution of Matthias Christian Schabel to Mathworks
Matthias Schabel (2020). 3D Shepp-Logan phantom
https://www.mathworks.com/matlabcentral/fileexchange/9416-3d-shepp-logan-phantom


2D phantoms are constructed on the grid of nxn voxels, 3D phantoms nxnxn and n is a parameter of the computation
Specification is however given on continuous grid of [-1,1] x [-1, 1] or [-1,1] x [-1, 1] x [-1, 1]
Continuous grids are mapped to voxel grids

In this file the following set of phantoms is predefined:

SheppLogan2D phantom
====================
Phantom from doi.org/10.1109/TNS.1974.6499235

	A			a			b			x0			y0			phi

  2.00000	  0.69000	  0.92000	  0.00000	  0.00000	  0.00000
 -0.98000	  0.66240	  0.87400	  0.00000	 -0.01840	  0.00000
 -0.02000	  0.11000	  0.31000	  0.22000	  0.00000	-18.00000
 -0.02000	  0.16000	  0.41000	 -0.22000	  0.00000	 18.00000
  0.01000	  0.21000	  0.25000	  0.00000	  0.35000	  0.00000
  0.01000	  0.04600	  0.04600	  0.00000	  0.10000	  0.00000
  0.01000	  0.04600	  0.04600	  0.00000	 -0.10000	  0.00000
  0.01000	  0.04600	  0.02300	 -0.08000	 -0.60500	  0.00000
  0.01000	  0.02300	  0.02300	  0.00000	 -0.60500	  0.00000
  0.01000	  0.02300	  0.04600	  0.06000	 -0.60500	  0.00000


Octave2D phantom
================
Taken from phantom.m file function shepp_logan of Octave package image-2.12.0 https://octave.sourceforge.io/image/

	A			a			b			x0			y0			phi

  1.00000	  0.69000	  0.92000	  0.00000	  0.00000	  0.00000
 -0.98000	  0.66240	  0.87400	  0.00000	 -0.01840	  0.00000
 -0.02000	  0.11000	  0.31000	  0.22000	  0.00000	-18.00000
 -0.02000	  0.16000	  0.41000	 -0.22000	  0.00000	 18.00000
  0.01000	  0.21000	  0.25000	  0.00000	  0.35000	  0.00000
  0.01000	  0.04600	  0.04600	  0.00000	  0.10000	  0.00000
  0.01000	  0.04600	  0.04600	  0.00000	 -0.10000	  0.00000
  0.01000	  0.04600	  0.02300	 -0.08000	 -0.60500	  0.00000
  0.01000	  0.02300	  0.02300	  0.00000	 -0.60600	  0.00000
  0.01000	  0.02300	  0.04600	  0.06000	 -0.60500	  0.00000

Toft2D phantom
==============
Also mod_shepp_logan image-2.12.0, taken from page 201, section B.3 of Peter Toft Ph.D. thesis https://petertoft.dk/PhD/PeterToft_PhD_thesis_5.pdf

	A			a			b			x0			y0			phi

  1.00000	  0.69000	  0.92000	  0.00000	  0.00000	  0.00000
 -0.80000	  0.66240	  0.87400	  0.00000	 -0.01840	  0.00000
 -0.20000	  0.11000	  0.31000	  0.22000	  0.00000	-18.00000
 -0.20000	  0.16000	  0.41000	 -0.22000	  0.00000	 18.00000
  0.10000	  0.21000	  0.25000	  0.00000	  0.35000	  0.00000
  0.10000	  0.04600	  0.04600	  0.00000	  0.10000	  0.00000
  0.10000	  0.04600	  0.04600	  0.00000	 -0.10000	  0.00000
  0.10000	  0.04600	  0.02300	 -0.08000	 -0.60500	  0.00000
  0.10000	  0.02300	  0.02300	  0.00000	 -0.60600	  0.00000
  0.10000	  0.02300	  0.04600	  0.06000	 -0.60500	  0.00000

ToftSchabel3D phantom
=================
Specification of 10 ellipsoids given in Matthias Schabel (2020). 3D Shepp-Logan phantom contribution inspirred by Toft2D
Matthias Schabel (2020). 3D Shepp-Logan phantom
https://www.mathworks.com/matlabcentral/fileexchange/9416-3d-shepp-logan-phantom

	A			a			b			c			x0			y0			z0			phi			theta		psi

  1.00000	  0.69000	  0.92000	  0.81000	  0.00000	  0.00000	  0.00000	  0.00000	  0.00000	  0.00000
 -0.80000	  0.66240	  0.87400	  0.78000	  0.00000	 -0.01840	  0.00000	  0.00000	  0.00000	  0.00000
 -0.20000	  0.11000	  0.31000	  0.22000	  0.22000	  0.00000	  0.00000	-18.00000	  0.00000	 10.00000
 -0.20000	  0.16000	  0.41000	  0.28000	 -0.22000	  0.00000	  0.00000	 18.00000	  0.00000	 10.00000
  0.10000	  0.21000	  0.25000	  0.41000	  0.00000	  0.35000	 -0.15000	  0.00000	  0.00000	  0.00000
  0.10000	  0.04600	  0.04600	  0.05000	  0.00000	  0.10000	  0.25000	  0.00000	  0.00000	  0.00000
  0.10000	  0.04600	  0.04600	  0.05000	  0.00000	 -0.10000	  0.25000	  0.00000	  0.00000	  0.00000
  0.10000	  0.04600	  0.02300	  0.05000	 -0.08000	 -0.60500	  0.00000	  0.00000	  0.00000	  0.00000
  0.10000	  0.02300	  0.02300	  0.20000	  0.00000	 -0.60600	  0.00000	  0.00000	  0.00000	  0.00000
  0.10000	  0.02300	  0.04600	  0.20000	  0.06000	 -0.60500	  0.00000	  0.00000	  0.00000	  0.00000



OctaveSchabel3D phantom
=================
Specification of 10 ellipsoids given in Matthias Schabel (2020). 3D Shepp-Logan phantom contribution inspirred by Octave2D
Matthias Schabel (2020). 3D Shepp-Logan phantom
https://www.mathworks.com/matlabcentral/fileexchange/9416-3d-shepp-logan-phantom

	 A			 a			 b			 c			 x0			 y0			 z0			 phi		 theta		 psi

  1.00000	  0.69000	  0.92000	  0.81000	  0.00000	  0.00000	  0.00000	  0.00000	  0.00000	  0.00000
 -0.98000	  0.66240	  0.87400	  0.78000	  0.00000	 -0.01840	  0.00000	  0.00000	  0.00000	  0.00000
 -0.02000	  0.11000	  0.31000	  0.22000	  0.22000	  0.00000	  0.00000	-18.00000	  0.00000	 10.00000
 -0.02000	  0.16000	  0.41000	  0.28000	 -0.22000	  0.00000	  0.00000	 18.00000	  0.00000	 10.00000
  0.01000	  0.21000	  0.25000	  0.41000	  0.00000	  0.35000	 -0.15000	  0.00000	  0.00000	  0.00000
  0.01000	  0.04600	  0.04600	  0.05000	  0.00000	  0.10000	  0.25000	  0.00000	  0.00000	  0.00000
  0.01000	  0.04600	  0.04600	  0.05000	  0.00000	 -0.10000	  0.25000	  0.00000	  0.00000	  0.00000
  0.01000	  0.04600	  0.02300	  0.05000	 -0.08000	 -0.60500	  0.00000	  0.00000	  0.00000	  0.00000
  0.01000	  0.02300	  0.02300	  0.20000	  0.00000	 -0.60600	  0.00000	  0.00000	  0.00000	  0.00000
  0.01000	  0.02300	  0.04600	  0.20000	  0.06000	 -0.60500	  0.00000	  0.00000	  0.00000	  0.00000

Kak3D phantom
=============
From ISBN 9780127745602 with some uncertainity

	 A			 a			 b			 c			 x0			 y0			 z0			 phi		 theta		 psi

  2.00000	  0.69000	  0.92000	  0.90000	  0.00000	  0.00000	  0.00000	  0.00000	  0.00000	  0.00000
 -0.98000	  0.66240	  0.87400	  0.88000	  0.00000	  0.00000	  0.00000	  0.00000	  0.00000	  0.00000
 -0.02000	  0.41000	  0.16000	  0.21000	 -0.22000	  0.00000	 -0.25000	 108.00000	  0.00000	  0.00000
 -0.02000	  0.31000	  0.11000	  0.22000	  0.22000	  0.00000	 -0.25000	 72.00000	  0.00000	  0.00000
  0.02000	  0.21000	  0.25000	  0.50000	  0.00000	  0.35000	 -0.25000	  0.00000	  0.00000	  0.00000
  0.02000	  0.04600	  0.04600	  0.04600	  0.00000	  0.10000	 -0.25000	  0.00000	  0.00000	  0.00000
  0.01000	  0.04600	  0.02300	  0.02000	 -0.08000	 -0.65000	 -0.25000	  0.00000	  0.00000	  0.00000
  0.01000	  0.04600	  0.02300	  0.02000	  0.06000	 -0.65000	 -0.25000	 90.00000	  0.00000	  0.00000
  0.02000	  0.05600	  0.04000	  0.10000	  0.06000	 -0.10500	  0.62500	 90.00000	  0.00000	  0.00000
 -0.02000	  0.05600	  0.05600	  0.10000	  0.00000	  0.10000	  0.62500	  0.00000	  0.00000	  0.00000


Yu3D phantom
============
From doi.org/10.1117/12.559300

	 A			 a			 b			 c			 x0			 y0			 z0			 phi		 theta		 psi

  1.00000	  0.69000	  0.92000	  0.90000	  0.00000	  0.00000	  0.00000	  0.00000	  0.00000	  0.00000
 -0.80000	  0.66240	  0.87400	  0.88000	  0.00000	  0.00000	  0.00000	  0.00000	  0.00000	  0.00000
 -0.20000	  0.41000	  0.16000	  0.21000	 -0.22000	  0.00000	 -0.25000	 108.00000	  0.00000	  0.00000
 -0.20000	  0.31000	  0.11000	  0.22000	  0.22000	  0.00000	 -0.25000	 72.00000	  0.00000	  0.00000
  0.20000	  0.21000	  0.25000	  0.50000	  0.00000	  0.35000	 -0.25000	  0.00000	  0.00000	  0.00000
  0.20000	  0.04600	  0.04600	  0.04600	  0.00000	  0.10000	 -0.25000	  0.00000	  0.00000	  0.00000
  0.10000	  0.04600	  0.02300	  0.02000	 -0.08000	 -0.65000	 -0.25000	  0.00000	  0.00000	  0.00000
  0.10000	  0.04600	  0.02300	  0.02000	  0.06000	 -0.65000	 -0.25000	 90.00000	  0.00000	  0.00000
  0.20000	  0.05600	  0.04000	  0.10000	  0.06000	 -0.10500	  0.62500	 90.00000	  0.00000	  0.00000
 -0.20000	  0.05600	  0.05600	  0.10000	  0.00000	  0.10000	  0.62500	  0.00000	  0.00000	  0.00000



Koay3D phantom
=============
From doi.org/ISBN 9780127745602

	 A			 a			 b			 c			 x0			 y0			 z0			 phi		 theta		 psi

  2.00000	  0.69000	  0.92000	  0.90000	  0.00000	  0.00000	  0.00000	  0.00000	  0.00000	  0.00000
 -0.80000	  0.66240	  0.87400	  0.88000	  0.00000	  0.00000	  0.00000	  0.00000	  0.00000	  0.00000
 -0.20000	  0.41000	  0.16000	  0.21000	 -0.22000	  0.00000	 -0.25000	 108.00000	  0.00000	  0.00000
 -0.20000	  0.31000	  0.11000	  0.22000	  0.22000	  0.00000	 -0.25000	 72.00000	  0.00000	  0.00000
  0.20000	  0.21000	  0.25000	  0.50000	  0.00000	  0.35000	 -0.25000	  0.00000	  0.00000	  0.00000
  0.20000	  0.04600	  0.04600	  0.04600	  0.00000	  0.10000	 -0.25000	  0.00000	  0.00000	  0.00000
  0.10000	  0.04600	  0.02300	  0.02000	 -0.08000	 -0.65000	 -0.25000	  0.00000	  0.00000	  0.00000
  0.10000	  0.04600	  0.02300	  0.02000	  0.06000	 -0.65000	 -0.25000	 90.00000	  0.00000	  0.00000
  0.20000	  0.05600	  0.04000	  0.10000	  0.06000	 -0.10500	  0.62500	 90.00000	  0.00000	  0.00000
 -0.20000	  0.05600	  0.05600	  0.10000	  0.00000	  0.10000	  0.62500	  0.00000	  0.00000	  0.00000

ToftSchabelKulvait3D phantom
===========================
ToftSchabel3D with all zero theta psi and z0 so that center slice is Toft2D

	 A			 a			 b			 c			 x0			 y0			 z0			 phi		 theta		 psi

  1.00000	  0.69000	  0.92000	  0.81000	  0.00000	  0.00000	  0.00000	  0.00000	  0.00000	  0.00000
 -0.80000	  0.66240	  0.87400	  0.78000	  0.00000	 -0.01840	  0.00000	  0.00000	  0.00000	  0.00000
 -0.20000	  0.11000	  0.31000	  0.22000	  0.22000	  0.00000	  0.00000	-18.00000	  0.00000	  0.00000
 -0.20000	  0.16000	  0.41000	  0.28000	 -0.22000	  0.00000	  0.00000	 18.00000	  0.00000	  0.00000
  0.10000	  0.21000	  0.25000	  0.41000	  0.00000	  0.35000	  0.00000	  0.00000	  0.00000	  0.00000
  0.10000	  0.04600	  0.04600	  0.05000	  0.00000	  0.10000	  0.00000	  0.00000	  0.00000	  0.00000
  0.10000	  0.04600	  0.04600	  0.05000	  0.00000	 -0.10000	  0.00000	  0.00000	  0.00000	  0.00000
  0.10000	  0.04600	  0.02300	  0.05000	 -0.08000	 -0.60500	  0.00000	  0.00000	  0.00000	  0.00000
  0.10000	  0.02300	  0.02300	  0.20000	  0.00000	 -0.60600	  0.00000	  0.00000	  0.00000	  0.00000
  0.10000	  0.02300	  0.04600	  0.20000	  0.06000	 -0.60500	  0.00000	  0.00000	  0.00000	  0.00000



The logic of rotations by theta and psi angles is taken from the code of Matthias Schabel. Other phantoms typically set theta=0 and psi=0 so there are no problems with different meanings of these angles in different authors.


Example use:
from denpy import PHANTOM
from denpy import DEN
from matplotlib import pyplot as plt

p = PHANTOM.phantom3d(n=256)
plt.imshow(p[:, :, 128], origin="lower")
plt.show()

#Exporting to DEN as FLOAT64
PHANTOM.storePhantom3DAsDEN("ToftSchabelKulvait3D_512.vol", phantom='ToftSchabelKulvait3D', n=512)


#Exporting to DEN as FLOAT32
from denpy import PHANTOM
from denpy import DEN
import numpy as np
p = PHANTOM.phantom3d(phantom='ToftSchabelKulvait3D', n=512)
p = p.astype(np.float32)
DEN.storeNdarrayAsDEN("ToftSchabelKulvait3D_512.vol", p, force=True)

Added code to create gauss shaped phantom.
"""
import numpy as np
import os
from collections import namedtuple
from denpy import DEN

Ellipsoid = namedtuple("Ellipsoid", "A a b c x0 y0 z0 phi theta psi")
Ellipse = namedtuple("Ellipse", "A a b x0 y0 phi")


def storePhantom3DAsDEN(fileName, phantom='ToftSchabelKulvait3D', n=512, force=False):
	if not force and os.path.exists(fileName):
		raise IOError('File already exists, no data written')
	phantom = phantom3d(phantom, n)
	DEN.storeNdarrayAsDEN(fileName, phantom, force=True)


def constructPhantom3D(e3d, n=512):
	p = np.zeros(n**3)
	rng = np.linspace(-1, 1, n)
	# The problem with the following is that when we later want to index
	# x[xind, yind, zind] so that x[xind,:,:] is constant as we apparently do
	# here then we have to use correct indexing moreover as in imagej axis y goes from top down, its better to flip it for better visualization
	x, y, z = np.meshgrid(rng, -rng, rng, indexing='ij')
	# Here we depend on internal alignment of numpy of flattened arrays and hope that flatten can be reversed by reshape
	# Note that for x[yind, xind, zind] indexing we would have to switch x and
	# y and flipping x, y, z in meshgrid function would lead to more natural
	# numpy alignment
	coord = np.vstack((x.flatten(), y.flatten(), z.flatten()))
	p = p.flatten()
	for e in e3d:
		phi = np.radians(e.phi)  # first Euler angle in radians
		theta = np.radians(e.theta)  # second Euler angle in radians
		psi = np.radians(e.psi)  # third Euler angle in radians
		c1 = np.cos(phi)
		s1 = np.sin(phi)
		c2 = np.cos(theta)
		s2 = np.sin(theta)
		c3 = np.cos(psi)
		s3 = np.sin(psi)
		# Euler rotation matrix
		# Here Matthias Schabel evidently has choosen so called Z1X2Z3 matrix as defined https://en.wikipedia.org/wiki/Euler_angles#Rotation_matrix
		# Implicit equation for ellipsoid is evaluated in rotated coordinates
		alpha = np.array([[c3 * c1 - c2 * s1 * s3, c3 * s1 + c2 * c1 * s3, s3 * s2],
		                  [-s3 * c1 - c2 * s1 * c3, -s3 * s1 + c2 * c1 * c3, c3 * s2], [s2 * s1, -s2 * c1, c2]])
		# rotated ellipsoid coordinates
		coordp = np.matmul(alpha, coord)
		ellipsoidCenter = np.array([[e.x0], [e.y0], [e.z0]])
		ellipsoidSquareHalfaxes = np.array(np.square([[e.a], [e.b], [e.c]]))
		ellipsoidImplicitResult = np.sum(np.divide(np.square(np.subtract(coordp, ellipsoidCenter)),
		                                           ellipsoidSquareHalfaxes),
		                                 axis=0)  # Sum along columns
		p[ellipsoidImplicitResult <= 1.0] += e.A
		"""
		idx = np.nonzero((coordp[0, :] - x0) ** 2.0 / asq +
						 (coordp[1, :] - y0) ** 2.0 / bsq +
						 (coordp[2, :] - z0) ** 2.0 / csq <= 1)[0]
		p[idx] = p[idx] + A
		"""
	p = p.reshape(n, n, n)
	# Here is the problem with default numpy indexing and data alignment, when
	# we index by [xind,yind,zind] but in numpy is expected to be indexed by [zind,yind,xind]
	# so that we swap it
	p = np.swapaxes(p, 0, 2)
	# Here we changed DEN.storeNdarrayAsDEN(fileName, phantom, force=True)indexing to [zind, yind, xind] as expected by DEN.storeNdarrayAsDEN that expects data in this format
	return p


def phantom3d(phantom='ToftSchabelKulvait3D', n=512):
	"""Three-dimensional Shepp-Logan like phantoms

Can be used to test 3D reconstruction algorithms.

Parameters
==========
	phantom: str
			One of {'modified-shepp-logan', 'shepp_logan', 'yu_ye_wang'},
			The type of phantom to draw.
	n : int, optional
			The grid size of the phantom

	Notes
	=====
	For any given voxel in the output image, the voxel's value is equal to the
	sum of the additive intensity values of all ellipsoids that the voxel is a
	part of.  If a voxel is not part of any ellipsoid, its value is 0.

	The additive intensity value A for an ellipsoid can be positive or
	negative;  if it is negative, the ellipsoid will be darker than the
	surrounding pixels.
	Note that, depending on the values of A, some voxels may have values
	outside the range [0,1].


Copyright 2005 Matthias Christian Schabel (matthias @ stanfordalumni . org)
University of Utah Department of Radiology
Utah Center for Advanced Imaging Research
729 Arapeen Drive
Salt Lake City, UT 84108-1218

This code is released under the Gnu Public License (GPL). For more information,
see : http://www.gnu.org/copyleft/gpl.html

	Source modified by Vojtech Kulvait
	"""
	if phantom == 'ToftSchabel3D':
		e3d = ToftSchabel3D()
	elif phantom == 'OctaveSchabel3D':
		e3d = OctaveSchabel3D()
	elif phantom == 'Kak3D':
		e3d = Kak3D()
	elif phantom == 'Yu3D':
		e3d = Yu3D()
	elif phantom == 'Koay3D':
		e3d = Koay3D()
	elif phantom == 'ToftSchabelKulvait3D':
		e3d = ToftSchabelKulvait3D()
	elif phantom == 'Chessboard3D':
		return Chessboard3D(n)
	else:
		raise TypeError('phantom type "%s" not recognized' % phantom)
	return constructPhantom3D(e3d, n)


def printEllipses(e2d, delimiter="\t"):
	print("A%sa%sb%sx0%sy0%sphi\n" % (delimiter, delimiter, delimiter, delimiter, delimiter),)
	for ellipse in e2d:
		print(delimiter.join("% 9.5f" % i for i in ellipse))


def printEllipsoids(e3d, delimiter="\t"):
	print(
	    "A%sa%sb%sc%sx0%sy0%sz0%sphi%stheta%spsi\n" %
	    (delimiter, delimiter, delimiter, delimiter, delimiter, delimiter, delimiter, delimiter, delimiter),)
	for ellipsoid in e3d:
		print(delimiter.join("% 9.5f" % i for i in ellipsoid))


def EllipseToEllipsoid(ellipse):
	"""Converts Ellipse to Ellipsoid where c=0, z0=0, theta=0, psi=0"""
	return Ellipsoid(ellipse.A, ellipse.a, ellipse.b, 0.0, ellipse.x0, ellipse.y0, 0.0, ellipse.phi, 0.0, 0.0)


def SheppLogan2D():
	"""
Shepp, L. A. & Logan, B. F.
The Fourier reconstruction of a head section
IEEE Transactions on Nuclear Science, Institute of Electrical and Electronics Engineers (IEEE), 1974, 21, 21-43
	"""
	e2d = []
	e2d.append(Ellipse(2.0, 0.69, 0.92, 0.0, 0.0, 0.0))
	e2d.append(Ellipse(-0.98, 0.6624, 0.874, 0.0, -0.0184, 0.0))
	e2d.append(Ellipse(-0.02, 0.11, 0.31, 0.22, 0.0, -18))
	e2d.append(Ellipse(-0.02, 0.16, 0.41, -0.22, 0.0, 18))
	e2d.append(Ellipse(0.01, 0.21, 0.25, 0.0, 0.35, 0.0))
	e2d.append(Ellipse(0.01, 0.046, 0.046, 0.0, 0.1, 0.0))
	e2d.append(Ellipse(0.01, 0.046, 0.046, 0.0, -0.1, 0.0))
	e2d.append(Ellipse(0.01, 0.046, 0.023, -0.08, -0.605, 0.0))
	e2d.append(Ellipse(0.01, 0.023, 0.023, 0.0, -0.605, 0.0))
	e2d.append(Ellipse(0.01, 0.023, 0.046, 0.06, -0.605, 0.0))
	return e2d


def Octave2D():
	"""Phantom differs from SheppLogan2D in the amplitude of the first ellipse being 1.0 instead of 2.0 and 9th ellipse y0 is -0.606 instead of -0.605"""
	# The Octave source code of image package states that:
	#
	# Note that the first element of this matrix, the gray value for the first
	# ellipse (human skull), has a value of 1.0 even though the paper gives it a
	# a value of 2.0 (see Table 1 on page 32 and Figure 1 on page 34). This
	# change is so that the **head** intensity values appear in the range [0 1]
	# rather than the range [1 2].
	#
	# **The problem with this**
	#
	# The background still need an intensity value which is going to be 0. This
	# means that we can't distinguish between the background and the ventricles
	# (ellipse "c" and "d" whose intensities are a + b + c and a + b + d, see
	# Figure 1) since they will have an intensity value of 0 (actually, because
	# of machine precision the ventricules will be almost 0). But if we didn't
	# made this change, the ** image** range would be [0 2] with all of the head
	# details compressed in half of the display range. Also, Matlab seems to be
	# doing the same.
	e2d = SheppLogan2D()
	e2d[0] = e2d[0]._replace(A=1.0)
	e2d[8] = e2d[8]._replace(y0=-0.606)
	return e2d


def Toft2D():
	"""Modified Shepp Logan from Octave"""
	e2d = SheppLogan2D()
	e2d[0] = e2d[0]._replace(A=1.0)
	e2d[1] = e2d[1]._replace(A=-0.8)
	e2d[2] = e2d[2]._replace(A=-0.2)
	e2d[3] = e2d[3]._replace(A=-0.2)
	e2d[4] = e2d[4]._replace(A=0.1)
	e2d[5] = e2d[5]._replace(A=0.1)
	e2d[6] = e2d[6]._replace(A=0.1)
	e2d[7] = e2d[7]._replace(A=0.1)
	e2d[8] = e2d[8]._replace(A=0.1, y0=-0.606)
	e2d[9] = e2d[9]._replace(A=0.1)
	return e2d


def ToftSchabel3D():
	"""
Specification of 10 ellipsoids given in Matthias Schabel (2020). 3D Shepp-Logan phantom contribution from modified_shepp_logan function
Matthias Schabel (2020). 3D Shepp-Logan phantom
https://www.mathworks.com/matlabcentral/fileexchange/9416-3d-shepp-logan-phantom
	"""
	e2d = Toft2D()
	e3d = [EllipseToEllipsoid(x) for x in e2d]
	e3d[0] = e3d[0]._replace(c=0.81)
	e3d[1] = e3d[1]._replace(c=0.78)
	e3d[2] = e3d[2]._replace(c=0.22, psi=10.0)
	e3d[3] = e3d[3]._replace(c=0.28, psi=10.0)
	e3d[4] = e3d[4]._replace(c=0.41, z0=-0.15)
	e3d[5] = e3d[5]._replace(c=0.05, z0=0.25)
	e3d[6] = e3d[6]._replace(c=0.05, z0=0.25)
	e3d[7] = e3d[7]._replace(c=0.05)
	e3d[8] = e3d[8]._replace(c=0.2)
	e3d[9] = e3d[9]._replace(c=0.2)
	return e3d


def OctaveSchabel3D():
	"""
Specification of 10 ellipsoids given in Matthias Schabel (2020). 3D Shepp-Logan phantom contribution from shepp_logan function
Matthias Schabel (2020). 3D Shepp-Logan phantom
https://www.mathworks.com/matlabcentral/fileexchange/9416-3d-shepp-logan-phantom
	"""
	e2d = Octave2D()
	e3d = [EllipseToEllipsoid(x) for x in e2d]
	e3d[0] = e3d[0]._replace(c=0.81)
	e3d[1] = e3d[1]._replace(c=0.78)
	e3d[2] = e3d[2]._replace(c=0.22, psi=10.0)
	e3d[3] = e3d[3]._replace(c=0.28, psi=10.0)
	e3d[4] = e3d[4]._replace(c=0.41, z0=-0.15)
	e3d[5] = e3d[5]._replace(c=0.05, z0=0.25)
	e3d[6] = e3d[6]._replace(c=0.05, z0=0.25)
	e3d[7] = e3d[7]._replace(c=0.05)
	e3d[8] = e3d[8]._replace(c=0.2)
	e3d[9] = e3d[9]._replace(c=0.2)
	return e3d


def Kak3D():
	"""
Taken from the book, Chapter 3, page 102
Kak, A.C. & Slaney, M.
Principles of computerized tomographic imaging
IEEE Press, 1988
ISBN 0879421983

but here apparently ellipses e and f are the same
corrected using the claim from
Koay, C. G.; Sarlls, J. E. & Özarslan, E.
Three-dimensional analytical magnetic resonance imaging phantom in the Fourier domain
Magnetic Resonance in Medicine, Wiley, 2007, 58, 430-436

so it should be the same as
Kak, A. & Young, T.
Handbook of pattern recognition and image processing
Academic Press, 1986
ISBN 9780127745602
	"""
	e3d = []
	e3d.append(Ellipsoid(2.0, 0.69, 0.92, 0.9, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0))
	e3d.append(Ellipsoid(-0.98, 0.6624, 0.874, 0.88, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0))
	e3d.append(Ellipsoid(-0.02, 0.41, 0.16, 0.21, -0.22, 0.0, -0.25, 108.0, 0.0, 0.0))
	e3d.append(Ellipsoid(-0.02, 0.31, 0.11, 0.22, 0.22, 0.0, -0.25, 72.0, 0.0, 0.0))
	e3d.append(Ellipsoid(0.02, 0.21, 0.25, 0.5, 0.0, 0.35, -0.25, 0.0, 0.0, 0.0))
	#Unsure about this
	e3d.append(Ellipsoid(0.02, 0.046, 0.046, 0.046, 0.0, 0.1, -0.25, 0.0, 0.0, 0.0))
	e3d.append(Ellipsoid(0.01, 0.046, 0.023, 0.02, -0.08, -0.65, -0.25, 0.0, 0.0, 0.0))
	e3d.append(Ellipsoid(0.01, 0.046, 0.023, 0.02, 0.06, -0.65, -0.25, 90.0, 0.0, 0.0))
	e3d.append(Ellipsoid(0.02, 0.056, 0.04, 0.1, 0.06, -0.105, 0.625, 90.0, 0.0, 0.0))
	e3d.append(Ellipsoid(-0.02, 0.056, 0.056, 0.1, 0.0, 0.1, 0.625, 0.0, 0.0, 0.0))
	return e3d


def Yu3D():
	"""
Yu, H.; Ye, Y. & Wang, G.
Bonse, U. (Ed.)
Katsevich-type algorithims for variable radius spiral cone-beam CT
Developments in X-Ray Tomography IV, SPIE, 2004
	"""
	e3d = Kak3D()
	e3d[0] = e3d[0]._replace(A=1.0)
	e3d[1] = e3d[1]._replace(A=-0.8)
	e3d[2] = e3d[2]._replace(A=-0.2)
	e3d[3] = e3d[3]._replace(A=-0.2)
	e3d[4] = e3d[4]._replace(A=0.2)
	e3d[5] = e3d[5]._replace(A=0.2)
	e3d[6] = e3d[6]._replace(A=0.1)
	e3d[7] = e3d[7]._replace(A=0.1)
	e3d[8] = e3d[8]._replace(A=0.2)
	e3d[9] = e3d[9]._replace(A=-0.2)
	return e3d


def Koay3D():
	"""
Koay, C. G.; Sarlls, J. E. & Özarslan, E.
Three-dimensional analytical magnetic resonance imaging phantom in the Fourier domain
Magnetic Resonance in Medicine, Wiley, 2007, 58, 430-436
	"""
	e3d = Kak3D()
	e3d[0] = e3d[0]._replace(A=2.0)
	e3d[1] = e3d[1]._replace(A=-0.8)
	e3d[2] = e3d[2]._replace(A=-0.2)
	e3d[3] = e3d[3]._replace(A=-0.2)
	e3d[4] = e3d[4]._replace(A=0.2)
	e3d[5] = e3d[5]._replace(A=0.2)
	e3d[6] = e3d[6]._replace(A=0.1)
	e3d[7] = e3d[7]._replace(A=0.1)
	e3d[8] = e3d[8]._replace(A=0.2)
	e3d[9] = e3d[9]._replace(A=-0.2)
	return e3d


def ToftSchabelKulvait3D():
	"""Same as ToftShabel3D but z0 and psi are all zero	"""
	e3d = ToftSchabel3D()
	e3d[2] = e3d[2]._replace(psi=0.0)
	e3d[3] = e3d[3]._replace(psi=0.0)
	e3d[4] = e3d[4]._replace(z0=0.0)
	e3d[5] = e3d[5]._replace(z0=0.0)
	e3d[6] = e3d[6]._replace(z0=0.0)
	return e3d


def Chessboard3D(n=512, block_count=8, i=None, j=None, k=None):
	"""
    Create a 3D chessboard phantom, with an option to isolate a single hypervoxel.

    Parameters:
    - n: int, size of the grid (n x n x n).
    - block_count: int, number of blocks along each dimension.
    - i: int or None, optional block index in the x-direction to fill (if specified).
    - j: int or None, optional block index in the y-direction to fill (if specified).
    - k: int or None, optional block index in the z-direction to fill (if specified).

    Returns:
    - grid: 3D numpy array of size (n, n, n) with the chessboard pattern or isolated block(s).
    """
	# Calculate the block size
	block_size = (n + block_count - 1) // block_count
	grid = np.zeros((n, n, n), dtype=np.float32)
	# Handle fixed hypervoxel indices if provided
	i_range = [i * block_size] if i is not None else range(0, n, block_size)
	j_range = [j * block_size] if j is not None else range(0, n, block_size)
	k_range = [k * block_size] if k is not None else range(0, n, block_size)
	# Loop over each block in the grid
	for i_start in i_range:
		for j_start in j_range:
			for k_start in k_range:
				# Determine the end indices for the current block, ensuring they don't exceed bounds
				i_end = min(i_start + block_size, n)
				j_end = min(j_start + block_size, n)
				k_end = min(k_start + block_size, n)
				# Determine whether to fill the block with 0 or 1 based on the parity of the sum of the indices
				if (i_start // block_size + j_start // block_size + k_start // block_size) % 2 == 0:
					# Fill the block with ones
					grid[k_start:k_end, j_start:j_end, i_start:i_end] = 1.0
	return grid


def construct2DGaussianDecay(centerx=256, centery=256, sigmax=5, sigmay=5, sizex=512, sizey=512, sizez=512):
	x = np.zeros(shape=(sizey, sizex), dtype=np.float32)
	factor = 1 / (2 * np.pi * (sigmax * sigmay))
	twosigmaxsquared = float(2 * sigmax * sigmax)
	twosigmaysquared = float(2 * sigmax * sigmax)
	for i in range(sizex):
		for j in range(sizey):
			x[j, i] = factor * np.exp(
			    (-float(i - centerx)**2 / twosigmaxsquared - float(j - centery)**2 / twosigmaysquared))
	return np.tile(x, (sizez, 1, 1))
