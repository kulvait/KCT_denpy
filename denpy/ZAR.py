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

def get_compressor(name, clevel=5, zarrv2=False, dtype=None, **codec_kwargs):
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
    codec_kwargs : dict
	    Additional keyword arguments to pass to the codec constructor.
        For example, ``bitspersample=12`` for codecs that support it.
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
			return Blosc(cname='zstd', clevel=int(clevel), shuffle=Blosc.BITSHUFFLE, typesize=itemsize)
		elif name == 'lz4' or name == 'blosc-lz4':
			return Blosc(cname='lz4', clevel=clevel, shuffle=Blosc.BITSHUFFLE, typesize=itemsize)
		elif name == 'gzip' or name == 'blosc-zlib':
			return GZip(level=clevel)
		elif name == 'blosc' or name == 'blosc-blosclz':
			return Blosc(cname='blosclz', clevel=clevel, shuffle=Blosc.BITSHUFFLE, typesize=itemsize)
		elif name == "avif":
			#clevel -1 AVIF_QUALITY_DEFAULT, 100 = AVIF_QUALITY_BEST = AVIF_QUALITY_LOSSLESS, 0 = AVIF_QUALITY_WORST
			from imagecodecs.numcodecs import Avif, register_codecs
			register_codecs()
			return Avif(level=clevel, **codec_kwargs)
		elif name == "jpegxr":
			from imagecodecs.numcodecs import Jpegxr, register_codecs
			register_codecs()
			return Jpegxr(level=clevel, **codec_kwargs)
		elif name == "jpegxl":
			from imagecodecs.numcodecs import Jpegxl, register_codecs
			register_codecs()
			if clevel == 0:
				jpegxl_codec = Jpegxl(lossless=True, **codec_kwargs)
			else:
				# Map clevel (typically 0-9) to JPEG XL distance parameter
				# clevel 1 = high quality, clevel 9 = low quality
				# distance: 0=lossless, 0.1-15=lossy (lower distance = higher quality)
				distance = max(0.1, (clevel - 1) * 1.5)  # scale clevel to distance
				jpegxl_codec = Jpegxl(lossless=False, distance=distance, effort=7, **codec_kwargs)
			return jpegxl_codec
		elif name == "jpeg2k":
			from imagecodecs.numcodecs import register_codecs, get_codec, Jpeg2k
			register_codecs()  # Ensure the codec is registered
			if clevel == 0:
				jp2_codec = get_codec({"id": Jpeg2k.codec_id, "reversible": True, **codec_kwargs})
			else:
				jp2_codec = get_codec({"id": Jpeg2k.codec_id, "reversible": False, "level": clevel, **codec_kwargs})
			return jp2_codec
		elif name == "htj2k":
			from imagecodecs.numcodecs import register_codecs, get_codec, Htj2k
			register_codecs()  # Ensure the codec is registered
			if clevel == 0:
				htj2k_codec = get_codec({"id": Htj2k.codec_id, "reversible": True, "level": clevel, **codec_kwargs})
			else:
				htj2k_codec = get_codec({"id": Htj2k.codec_id, "reversible": False, **codec_kwargs})
			return htj2k_codec
		elif name == "sz3":
			import imagecodecs
			if not imagecodecs.SZ3.available:
				raise RuntimeError("SZ3 codec is not available in the current Python imagecodecs package.")
			from imagecodecs.numcodecs import register_codecs, get_codec, Sz3
			register_codecs()  # Ensure the codec is registered
			return get_codec({"id": Sz3.codec_id, "mode": "abs", "abs": clevel})  # Use clevel as absolute error for lossy compression
		elif name == "zfp":
			import imagecodecs
			if not imagecodecs.ZFP.available:
				raise RuntimeError("ZFP codec is not available in the current Python imagecodecs package.")
			from imagecodecs.numcodecs import register_codecs, get_codec, Zfp
			register_codecs()  # Ensure the codec is registered
			if clevel == 0:
				# For lossless compression, use mode="lossless"
				return get_codec({"id": Zfp.codec_id, "mode": imagecodecs.ZFP.MODE.REVERSIBLE, "numthreads": os.cpu_count()})  # Use clevel as absolute error for lossy compression
			else:
				return get_codec({"id": Zfp.codec_id, "mode": imagecodecs.ZFP.MODE.FIXED_PRECISION, "level": clevel, "numthreads": os.cpu_count()})  # Use clevel as absolute error for lossy compression
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
		elif name == "lz4":
			codecs_chain.append(codecs.LZ4Codec(level=clevel))
		elif name == "gzip":
			codecs_chain.append(codecs.GzipCodec(level=clevel))
		elif name == "avif":
			from imagecodecs.zarr import Avif, register_codecs
			register_codecs()
			# clevel -1 AVIF_QUALITY_DEFAULT, 100 = AVIF_QUALITY_BEST = AVIF_QUALITY_LOSSLESS, 0 = AVIF_QUALITY_WORST
			# Note current implementation for momochrome always force AVIF_QUALITY_LOSSLESS, so clevel is ignored for monochrome images
			codecs_chain.append(Avif(level=clevel, **codec_kwargs))
		elif name == "jpegxr":
			from imagecodecs.zarr import Jpegxr, register_codecs
			register_codecs()
			Jpegxr_codec = Jpegxr(level=clevel, **codec_kwargs)
			codecs_chain.append(Jpegxr_codec)
		elif name == "jpegxl":
			from imagecodecs.zarr import Jpegxl
			# JPEG XL offers both lossless and lossy compression
			#  # -inf-100: quality; > 100: lossless
			if clevel == 0:
				jpegxl_codec = Jpegxl(lossless=True, **codec_kwargs)
			else:
				distance = max(0.1, (clevel - 1) * 1.5)  # scale clevel to distance
				jpegxl_codec = Jpegxl(lossless=False, distance=distance, effort=7, **codec_kwargs)
			codecs_chain.append(jpegxl_codec)
		elif name == "jpeg2k":
			from imagecodecs.zarr import Jpeg2k
			#"Zarr v3: Using JPEG 2000 codec with clevel={clevel}. Note: For lossless compression, use clevel=0."
			# quality, psnr, level < 1 or > 1000 map to quality=0
			if clevel == 0:
				jp2_codec = Jpeg2k(reversible=True, **codec_kwargs)
			else:
				jp2_codec = Jpeg2k(reversible=False, level=clevel, **codec_kwargs)
			codecs_chain.append(jp2_codec)
			# Try level 5, for lossy implementation, use reversible=False
		elif name == "htj2k":
			from imagecodecs.zarr import Htj2k
			# Note that level maps to qstep [0.0-1.0) or qfactor [1-100] ... qfactor 99 too low quatlity, so use qstep=clevel=(0,0.5) for lossy compression lower is better quality, 0.0=lossless, 0.5=highly lossy
			if clevel == 0:
				htj2k_codec = Htj2k(reversible=True, **codec_kwargs)
			else:
				htj2k_codec = Htj2k(reversible=False, level=clevel, **codec_kwargs)
			codecs_chain.append(htj2k_codec)
		elif name == "sz3":
			import imagecodecs
			if not imagecodecs.SZ3.available:
				raise RuntimeError("SZ3 codec is not available in the current Python imagecodecs package.")
			from imagecodecs.zarr import Sz3
			sz3_codec = Sz3(mode="abs", abs=clevel)
			codecs_chain.append(sz3_codec)
		elif name == "zfp":
			import imagecodecs
			if not imagecodecs.ZFP.available:
				raise RuntimeError("ZFP codec is not available in the current Python imagecodecs package.")
			from imagecodecs.zarr import Zfp
			if clevel == 0:
				# For lossless compression, use mode="lossless"
				zfp_codec = Zfp(mode=imagecodecs.ZFP.MODE.REVERSIBLE, numthreads=os.cpu_count())
			else:
				zfp_codec = Zfp(mode=imagecodecs.ZFP.MODE.FIXED_PRECISION, level=clevel, numthreads=os.cpu_count())
			codecs_chain.append(zfp_codec)
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
		elif name == "blosc-zstd" or name=="zstd":
			codecs_chain.append(
				codecs.BloscCodec(
					cname=codecs.BloscCname.zstd,
					clevel=int(clevel),
					shuffle="shuffle",
					typesize=itemsize,
				)
			)
		else:
			raise ValueError(f"Unknown compressor type '{name}' for Zarr v3")
		return codecs_chain
