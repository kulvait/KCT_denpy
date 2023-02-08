#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
(c) 2019

DEN and DICOM IO manipulation

@author: Vojtech Kulvait
@author: Enno Schieferecke

"""

import glob
import os
import pydicom

# https://stackoverflow.com/questions/46304306/how-to-generate-unique-dicom-uid
def generateRandomUID(prefix = None):
	import uuid
	id128 = uuid.uuid4().int
	if prefix is None:
		return  "2.25.%d"%(id128)
	else:
		prefix = str(prefix)
		prefixlen = len(prefix)
		suffixlen = 64 - prefixlen
		suffix = "%d"%id128
		suffix = suffix[:suffixlen]
		uid = "%s%s"%(prefix, suffix)[:64]
		return uid

# In dicom the time is encoded by type TM
# A string of characters of the format HHMMSS.FFFFFF; where HH contains
# hours (range "00" - "23"), MM contains minutes (range "00" - "59"), SS
# contains seconds (range "00" - "60"), and FFFFFF contains a fractional
# part of a second as small as 1 millionth of a second (range "000000" -
# "999999"). A 24-hour clock is used. Midnight shall be represented by
# only "0000" since "2400" would violate the hour range. The string may be
# padded with trailing spaces. Leading and embedded spaces are not
# allowed.
def timeToSeconds(timestring):
	import datetime
	timestring = timestring.strip();
	if "." in timestring:
		dt = datetime.datetime.strptime(timestring, "%H%M%S.%f")
	else:
		dt = datetime.datetime.strptime(timestring, "%H%M%S")
	ddt = dt - datetime.datetime(1900, 1, 1)
	return(ddt.total_seconds())

# datestring contains DA and timestring DT
# http://dicom.nema.org/medical/dicom/current/output/chtml/part05/sect_6.2.html


def dateAndTimeToSeconds(datestring, timestring):
	datestring = datestring.strip()
	timestring = timestring.strip()
	return(datetimeToSeconds(datestring + timestring))

# datestimetring contains DT from http://dicom.nema.org/medical/dicom/current/output/chtml/part05/sect_6.2.html
# A concatenated date-time character string in the format:
# YYYYMMDDHHMMSS.FFFFFF&ZZXX
# The components of this string, from left to right, are YYYY = Year, MM = Month, DD = Day, HH = Hour (range "00" - "23"), MM = Minute (range "00" - "59"), SS = Second (range "00" - "60").
# FFFFFF = Fractional Second contains a fractional part of a second as small as 1 millionth of a second (range "000000" - "999999").
#&ZZXX is an optional suffix for offset from Coordinated Universal Time (UTC), where & = "+" or "-", and ZZ = Hours and XX = Minutes of offset.
def datetimeToSeconds(datetimestring):
	datetimestring=datetimestring.strip()
	if datetimestring.find("&") != -1:	# Remove all after
		datetimestring = datetimestring[:datetimestring.find("&")]
	import datetime
	if "." in datetimestring:
		dt = datetime.datetime.strptime(datetimestring, "%Y%m%d%H%M%S.%f")
	else:
		dt = datetime.datetime.strptime(datetimestring, "%Y%m%d%H%M%S")
	ddt = dt - datetime.datetime(1900, 1, 1)
	return(ddt.total_seconds())


def getDicomFiles(d, suffixes=["IMA"]):
	dicomFiles = []
	for s in suffixes:
		dicomFiles.extend(glob.glob(os.path.join(d, "*.%s" % s)))
	dicomFiles.sort()
	return dicomFiles

# Get list of pydicom objects for each file with given suffix in the directory


def getDicoms(d, suffixes=["IMA"]):
	dicomFiles = getDicomFiles(d, suffixes)
	dicoms = [pydicom.read_file(x) for x in dicomFiles]
	return dicoms


def getPtsIds(dicoms):
	patientIDs = list(
		set([x.PatientID for x in dicoms if "PatientID" in x.dir("PatientID")]))
	return patientIDs


def getStudyDates(dicoms):
	patientIDs = list(
		set([x.PatientID for x in dicoms if "StudyDate" in x.dir("StudyDate")]))
	return patientIDs


def getSeriesDescriptions(dicoms):
	IDs = list(
		set([x.SeriesDescription for x in dicoms if "SeriesDescription" in x.dir("SeriesDescription")]))
	return IDs


def getInstanceUIDs(dicoms):
	IDs = list(
		set([x.SeriesInstanceUID for x in dicoms if "SeriesDescription" in x.dir("SeriesDescription")]))
	return IDs
