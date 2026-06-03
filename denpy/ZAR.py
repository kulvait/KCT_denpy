#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created 04 2026

Zarr format manipulation and helper functions.

@author: Vojtěch Kulvait
@license: GNU General Public License v3.0

"""
import zarr
import os
import numpy as np
import logging


# Create a logger specific to this module
log = logging.getLogger(__name__)
log.setLevel(logging.INFO) # Set the logging level to INFO
# Create a console handler and set its level to INFO
ch = logging.StreamHandler()
ch.setLevel(logging.INFO)
# Create a formatter and set it for the handler
formatter = logging.Formatter('%(asctime)s - %(name)s:%(lineno)d - %(levelname)s : %(message)s', datefmt='%d.%m.%Y %H:%M:%S')
ch.setFormatter(formatter)
# Add the handler to the logger
log.addHandler(ch)
log.propagate = False  # Prevent log messages from being propagated to the root logger

#np.arrays are by default in row-major order and they are indexed as follows
#array2d.shape = (dimy, dimx) = (axis0, axis1)
#array3d.shape=(dimz, dimy, dimx)= (axis0, axis1, axis2)

def get_compressor(name, clevel=5, zarrv2=False, dtype=None):
	"""
	Return a zarr-compatible compressor/codec based on name and Zarr format version.

	Parameters
	----------
	name : str
		Compression name (e.g., 'none', 'zstd', 'blosc-zstd', 'lz4', 'gzip', ...).
	clevel : int
		Compression level (meaning depends on the codec; Zstd/Blosc: 0..9 typical).
	zarrv2 : bool
		If False, return a Zarr v3 codec *pipeline* (list) suitable for `codecs=...`.
		If True, return a single compressor object (e.g., for Zarr v2 `compressor=`).
		Default is False (Zarr v3).
	dtype : Optional[Union[np.dtype, type, str]]
		Array dtype (e.g., np.uint16, 'uint16', np.dtype('uint16')). Used to set
		Blosc `typesize` (bytes per element). If None, defaults to itemsize=1.
		Important for shuffle codecs like Blosc, which require a typesize to function correctly.
	"""
	# Derive typesize from outtype if provided
	itemsize = 1  # Default typesize for codecs that require it (e.g., Blosc)
	if dtype is not None:
		if isinstance(dtype, np.dtype):
			typesize = dtype.itemsize
		else:
			try:
				typesize = np.dtype(dtype).itemsize
			except TypeError:
				log.warning(f"Invalid dtype provided: {dtype}. Defaulting to typesize=1.")
	if zarrv2:
		# ---- Zarr v2 codec  ----
		from numcodecs import Blosc, GZip as NcGZip
		# Old style compressors (zarr v2 compatible)
		if name == 'none':
			return None
		elif name == 'zstd' or name == 'blosc-zstd':
			return Blosc(cname='zstd', clevel=clevel, shuffle=Blosc.BITSHUFFLE, typesize=itemsize)
		elif name == 'lz4' or name == 'blosc-lz4':
			return Blosc(cname='lz4', clevel=clevel, shuffle=Blosc.BITSHUFFLE, typesize=itemsize)
		elif name == 'gzip' or name == 'blosc-zlib':
			return GZip(level=clevel)
		elif name == 'blosc' or name == 'blosc-blosclz':
			return Blosc(cname='blosclz', clevel=clevel, shuffle=Blosc.BITSHUFFLE, typesize=itemsize)
		elif name == "avif":
			from imagecodecs.numcodecs import Avif, register_codecs
			register_codecs()
			return Avif(bitspersample=12, numthreads=os.cpu_count())  # Return the codec instance directly for Zarr v2
		elif name == "jpegxr":
			from imagecodecs.numcodecs import Jpegxr
			return Jpegxr()  # Return the codec instance directly for Zarr v2
		elif name == "jpeg2k":
			from imagecodecs.numcodecs import register_codecs, get_codec, Jpeg2k
			register_codecs()  # Ensure the codec is registered
			if clevel == 0:
				jp2_codec = get_codec({"id": Jpeg2k.codec_id, "bitspersample": 12, "reversible": True, "colorspace": "GRAY", "mct": False, "numthreads": os.cpu_count()})
			else:
				#jp2_codec = get_codec({"id": Jpeg2k.codec_id, "bitspersample": 12, "reversible": False, "colorspace": "GRAY", "mct": False, "level": clevel})
				jp2_codec = get_codec({"id": Jpeg2k.codec_id, "reversible": False, "colorspace": "GRAY", "mct": False, "level": clevel, "numthreads": os.cpu_count()})
			return jp2_codec
		else:
			raise ValueError(f"Unknown compression type: {name}")
	else:
		# ---- Zarr v3 codecs (lazy import for safety) ----
		try:
			import zarr.codecs as codecs
		except ImportError:
			raise ImportError(
				"Zarr v3 codec system not available in this version of zarr. "
				"Please upgrade to zarr>=2.18.0."
			)
		# Map names to codecs
		codecs_chain = []
		if name == 'none':
			print("No compression selected for Zarr v3, returning empty codec chain.")
		elif name == "zstd":
			codecs_chain.append(codecs.ZstdCodec(level=clevel))
		elif name == "lz4":
			codecs_chain.append(codecs.LZ4Codec(level=clevel))
		elif name == "gzip":
			codecs_chain.append(codecs.GzipCodec(level=clevel))
		elif name == "avif":
			from imagecodecs.numcodecs import register_codecs, Avif
			register_codecs(Avif)  # Ensure the codec is registered
			codecs_chain.append(zarr.get_codec({"id": Avif.codec_id}))
		elif name == "jpegxr":
			from imagecodecs.numcodecs import Jpegxr
			register_codecs()
			codecs_chain.append(zarr.get_codec({"id": Jpegxr.codec_id}))
		elif name == "jpeg2k":
			#For Zarr v3 there is still no functional JPEG 2000 codec, but we can use the imagecodecs implementation as a custom codec in the chain
			from imagecodecs.zarr import Jpeg2k
			#"Zarr v3: Using JPEG 2000 codec with clevel={clevel}. Note: For lossless compression, use clevel=0."
			if clevel == 0:
				jp2_codec = Jpeg2k(reversible=True, colorspace="GRAY", mct=False, bitspersample=12, numthreads=os.cpu_count())
			else:
				jp2_codec = Jpeg2k(reversible=False, colorspace="GRAY", mct=False, bitspersample=12, level=clevel, numthreads=os.cpu_count())
			codecs_chain.append(jp2_codec)
			# Try level 5, for lossy implementation, use reversible=False
		elif name == "blosc" or name == "blosc-blosclz":
			codecs_chain.append(
				codecs.BloscCodec(
					cname=codecs.BloscCname.blosclz,
					clevel=clevel,
					shuffle="shuffle",
					typesize=itemsize,
				)
			)
		elif name == "blosc-lz4":
			codecs_chain.append(
				codecs.BloscCodec(
					cname=codecs.BloscCname.lz4,
					clevel=clevel,
					shuffle="shuffle",
					typesize=itemsize,
				)
			)
		elif name == "blosc-lz4hc":
			codecs_chain.append(
				codecs.BloscCodec(
					cname=codecs.BloscCname.lz4hc,
					clevel=clevel,
					shuffle="shuffle",
					typesize=itemsize,
				)
			)
		elif name == "blosc-snappy":
			codecs_chain.append(
				codecs.BloscCodec(
					cname=codecs.BloscCname.snappy,
					clevel=clevel,
					shuffle="shuffle",
					typesize=itemsize,
				)
			)
		elif name == "blosc-zlib":
			codecs_chain.append(
				codecs.BloscCodec(
					cname=codecs.BloscCname.zlib,
					clevel=clevel,
					shuffle="shuffle",
					typesize=itemsize,
				)
			)
		elif name == "blosc-zstd":
			codecs_chain.append(
				codecs.BloscCodec(
					cname=codecs.BloscCname.zstd,
					clevel=clevel,
					shuffle="shuffle",
					typesize=itemsize,
				)
			)
		else:
			raise ValueError(f"Unknown compressor type '{name}' for Zarr v3")
		return codecs_chain
