#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created 04 2019

Functions for ECC based image registration of DEN files

@author: Vojtech Kulvait

"""
import numpy as np
import os

# Estimate rigid transformation
# Most naive approach


def getRigidTransformation(frame1, frame2):
    minvalue = min(np.percentile(frame1, 10), np.percentile(frame2, 10))
    maxvalue = min(np.percentile(frame1, 90), np.percentile(frame2, 90))
    frame1 = (frame1 - minvalue) / (maxvalue - minvalue)
    frame2 = (frame2 - minvalue) / (maxvalue - minvalue)
    frame1 = np.maximum(0, np.minimum(1, frame1))
    frame2 = np.maximum(0, np.minimum(1, frame2))
    warp_matrix = cv2.estimateRigidTransform(
        (frame1 * 255 + 0.5).astype('uint8'), (frame2 * 255 + 0.5).astype('uint8'), False)
    return warp_matrix

# Estimate of corelation coefficient
# Gets result after first step
# There is no function for direct computation


def estimateCC(frame1, frame2):
    warp_mode = cv2.MOTION_EUCLIDEAN
    warp_matrix = np.eye(2, 3, dtype=np.float32)
    max_iterations = 1
    termination_eps = 1e-15
    criteria = (cv2.TERM_CRITERIA_EPS |
                cv2.TERM_CRITERIA_COUNT, max_iterations, termination_eps)
    (cc, warp_matrix) = cv2.findTransformECC(
        frame1, frame2, warp_matrix, warp_mode, criteria)
    return (cc, warp_matrix)

# Get optimal transformation


def getTransformECC(frame1, frame2, warp_mode=None, max_iterations=None, termination_eps=None, init_warp_matrix=None):
    if warp_mode is None:
        warp_mode = cv2.MOTION_EUCLIDEAN
    if max_iterations is None:
        max_iterations = 5000  # Default of findTransformECC is 50
    if termination_eps is None:
        termination_eps = 1e-10
    # Default of findTransformECC is 0.001 if init_warp_matrix is None
    if warp_mode == cv2.MOTION_HOMOGRAPHY:
        init_warp_matrix = np.eye(3, 3, dtype=np.float32)
    else:
        init_warp_matrix = np.eye(2, 3, dtype=np.float32)
    criteria = (cv2.TERM_CRITERIA_EPS |
                cv2.TERM_CRITERIA_COUNT, max_iterations, termination_eps)
    (cc, warp_matrix) = cv2.findTransformECC(frame1,
                                             frame2, init_warp_matrix, warp_mode, criteria)
    return (cc, warp_matrix)

# Calculate the x and y gradients using Sobel operator
# Then using Euclidean norm of both


def getGradient(frame):
    grad_x = cv2.Sobel(frame, cv2.CV_32F, 1, 0, ksize=3)
    grad_y = cv2.Sobel(frame, cv2.CV_32F, 0, 1, ksize=3)
    return np.sqrt(np.add(np.square(grad_x), np.square(grad_y)))

# Calculate Laplacian


def getLaplacian(frame):
    return (cv2.Laplacian(frame, cv2.CV_32F, ksize=23))

# Apply the same algorithm on gradient data


def gradientTransformECC(frame1, frame2, warp_mode=None, max_iterations=None, termination_eps=None, init_warp_matrix=None):
    g1 = getGradient(frame1)
    g2 = getGradient(frame2)
    return getTransformECC(g1, g2, warp_mode, max_iterations, termination_eps, init_warp_matrix)

# Apply the same algorithm on laplacian data


def laplacianTransformECC(frame1, frame2, warp_mode=None, max_iterations=None, termination_eps=None, init_warp_matrix=None):
    l1 = getLaplacian(frame1)
    l2 = getLaplacian(frame2)
    return getTransformECC(l1, l2, warp_mode, max_iterations, termination_eps, init_warp_matrix)

# For solving the problem that we have warp_matrix in the coordinates(j, i)
# here(x, y) = (col, row)
# We want the warp_matrix to operate on the data relative to(j ',i')
#(i, j) =(i '-offset_row, j' - offset_col)
# We have that(Ax) ' = Ax' - Aoffset + offset
# Here offset is(offset_row, offset_col, 0)
# So the operator that acts on(i ',j') space has modified third column.


def shiftWarpMatrix(warp_matrix, offset_row, offset_col):
    output_matrix = np.copy(warp_matrix)
    offsetVector = np.array([offset_col, offset_row, 0], dtype=np.float32)
    thirdRowAdd = - \
        np.matmul(warp_matrix, offsetVector) + \
        offsetVector[0:warp_matrix.shape[0]]
    output_matrix[:, 2] = output_matrix[:, 2] + thirdRowAdd
    return output_matrix

# Perform the 'backtransformation' of the frame to the template space
# It seems that this function takes(x, y) indexing(cols, rows)
# Interpolation mode see
#
# https://docs.opencv.org/3.1.0/da/d54/group__imgproc__transform.html#gga5bb5a1fea74ea38e1a5445ca803ff121aa5521d8e080972c762467c45f3b70e6c


def warpFrame(frame, warp_matrix, warp_mode, interpolation_mode=None):
    if interpolation_mode is None:
        interpolation_mode = cv2.INTER_LANCZOS4
        (rows, cols) = frame.shape
    if warp_mode == cv2.MOTION_HOMOGRAPHY:
        frameAligned = cv2.warpPerspective(frame, warp_matrix, (
            cols, rows), flags=interpolation_mode + cv2.WARP_INVERSE_MAP)
    else:
        frameAligned = cv2.warpAffine(frame, warp_matrix, (
            cols, rows), flags=interpolation_mode + cv2.WARP_INVERSE_MAP)
    return frameAligned
