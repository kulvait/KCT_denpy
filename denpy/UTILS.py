#/usr/bin/env python
#-*- coding: utf-8 -*-

#Author: Vojtech Kulvait
#Date: 2018-2024
#License: GPL3

#Utils to easily process files in the pipeline
import os


#Reads params file and possibly updates the dictionary params with its contents
def readParamsFile(f, params=None):
	if params is None:
		params = {}
	try:
		with open(f) as p:
			for l in p:
				l = l.strip()
				if not l.startswith("#") and "=" in l:
					name, var = l.partition("=")[::2]
					params[name.strip()] = var.strip()
		return (params)
	except IOError:
		print("Could not read params file %s" % f)
		return (params)


def writeParamsFile(params, f):
	with open(f, "w") as p:
		p.write("#!/bin/bash")
		p.write("\n")
		kx = list(params.keys())
		kx.sort()
		for k in kx:
			p.write("%s=%s" % (str(k), str(params[k])))
			p.write("\n")
		p.close()
