#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created 11 2026

Convert stacks of TIFF files to zarr or DEN as described by a JSON specification
(denpy docs/TIFFSTACK_SPEC.md), or build such a specification from a list of files.

	convertTiffStack scan.json [-j 16] [--report report.json] [--dry-run] [--check-headers]
	convertTiffStack --files img_*.tif --zarr img.zar [--compression htj2k --clevel 0] [--write-spec s.json]
	convertTiffStack --files img_*.tif --den img.den [--dtype float32]

@author: Vojtěch Kulvait
@license: GNU General Public License v3.0
"""
import argparse
import json
import os
import re
import sys
import time

from denpy import TIFF


def natural_key(text):
	"""img2 before img10."""
	return [int(t) if t.isdigit() else t for t in re.split(r"(\d+)", text)]


def build_parser():
	parser = argparse.ArgumentParser(
		prog="convertTiffStack",
		description="Convert stacks of TIFF files to zarr or DEN as described by a JSON specification "
		"(denpy docs/TIFFSTACK_SPEC.md), or convert a list of files to one array.")
	parser.add_argument("spec", nargs="?", help="JSON specification (not with --files)")
	parser.add_argument("-j", "--threads", type=int, default=-1,
						help="Worker threads. [defaults to -1: a quarter of the CPUs available, at least %d, at most %d]"
						% (TIFF.DEFAULT_MIN_THREADS, TIFF.DEFAULT_MAX_THREADS))
	parser.add_argument("--report", help="Write the report (per array: shape, dtype, null frames, bytes, time) as JSON.")
	parser.add_argument("--dry-run", action="store_true", help="Validate and print what would be written, write nothing.")
	parser.add_argument("--check-headers", action="store_true",
						help="Validation reads the header of every TIFF and compares shape and dtype with the spec.")
	parser.add_argument("--max-memory", type=float, default=None,
						help="Limit in GiB for reductions holding a whole stack (median). [defaults to a quarter of the RAM]")
	parser.add_argument("--extend", action="append", default=[], metavar="EXTENSION.json",
						help="Extend the spec (virtualstacks, arrays, defaults; see TIFF.extend_spec); can be repeated.")
	parser.add_argument("--verbose", action="store_true", help="Print every written block.")
	files = parser.add_argument_group("one array from a list of files (instead of a spec)")
	files.add_argument("--files", nargs="+", help="TIFF files, one frame each, in this order (see --sort).")
	target = files.add_mutually_exclusive_group()
	target.add_argument("--zarr", help="Output zarr store, the root is the array.")
	target.add_argument("--den", help="Output DEN file.")
	files.add_argument("--store", choices=("zip", "directory"), default="zip", help="Zarr store. [defaults to zip]")
	files.add_argument("--zarr-format", type=int, choices=(2, 3), default=3, help="[defaults to 3]")
	files.add_argument("--compression", default="none", help="Zarr codec name as in denpy.ZAR.get_compressor. [defaults to none]")
	files.add_argument("--clevel", type=float, default=None, help="Codec level (0 = lossless for htj2k, jpeg2k, jpegxl).")
	files.add_argument("--dtype", help="Output dtype. [defaults to the dtype of the files]")
	files.add_argument("--clip", action="store_true", help="Clip values outside the range of an integer --dtype.")
	files.add_argument("--sort", choices=("none", "natural"), default="none",
					   help="Order of --files: as given (shell glob order) or natural (img2 before img10). [defaults to none]")
	files.add_argument("--mode", choices=("w-", "w"), default="w-", help="w-: fail if the output exists, w: replace it. [defaults to w-]")
	files.add_argument("--write-spec", help="Save the specification (built from --files, or the spec with its --extend extensions).")
	return parser


def spec_from_arguments(arg, parser):
	if arg.spec:
		parser.error("give a spec or --files, not both")
	if not (arg.zarr or arg.den):
		parser.error("--files needs --zarr or --den")
	files = sorted(arg.files, key=natural_key) if arg.sort == "natural" else list(arg.files)
	if arg.zarr:
		output = {"format": "zarr", "store": arg.store, "path": arg.zarr, "zarr_format": arg.zarr_format, "mode": arg.mode}
	else:
		output = {"format": "den", "path": arg.den, "mode": arg.mode}
	array = {}
	if arg.dtype:
		array["dtype"] = arg.dtype
	if arg.clip:
		array["clip"] = True
	if arg.zarr and arg.compression != "none":
		codec = {"name": arg.compression}
		if arg.clevel is not None:
			codec["clevel"] = int(arg.clevel) if float(arg.clevel).is_integer() else arg.clevel
		array["zarr"] = {"codec": codec}
	return TIFF.spec_from_files(files, output, base_dir=os.getcwd(), **array)


def main(argv=None):
	parser = build_parser()
	argv = sys.argv[1:] if argv is None else list(argv)
	arg = parser.parse_args(argv or ["--help"])
	start = time.time()
	print("START convertTiffStack " + " ".join(argv), flush=True)
	try:
		if arg.files:
			spec = spec_from_arguments(arg, parser)
		elif arg.spec:
			spec, base_dir = TIFF.load_spec(arg.spec)
			spec["base_dir"] = base_dir  # resolved against the spec file, as when it is read by validate
		else:
			parser.error("give a spec or --files")
		for extension in arg.extend:
			spec = TIFF.extend_spec(spec, extension)
		if arg.write_spec:
			with open(arg.write_spec, "w") as f:
				json.dump(spec, f, indent=1)
				f.write("\n")
		max_memory = None if arg.max_memory is None else int(arg.max_memory * 2**30)
		plan = TIFF.validate(spec, check_headers=arg.check_headers, max_memory=max_memory)
		print(plan.describe(), flush=True)
		if arg.dry_run:
			return 0
		report = TIFF.execute(plan, threads=arg.threads, verbose=arg.verbose, log=lambda m: print(m, flush=True))
	except (TIFF.SpecError, TIFF.FrameError) as e:
		print(f"ERROR: {e}", file=sys.stderr)
		return 1
	for name, a in report["arrays"].items():
		print(f"SIZE {name or 'array'}: {a['stored_bytes'] / 2**30:.3f} GiB stored, "
			  f"{a['ratio_to_raw']:.3f} of the raw array, TIFF input {a['input_bytes'] / 2**30:.3f} GiB", flush=True)
	if arg.report:
		with open(arg.report, "w") as f:
			json.dump(report, f, indent=1)
			f.write("\n")
	print(f"END convertTiffStack {time.time() - start:.1f}s", flush=True)
	return 0


if __name__ == "__main__":
	sys.exit(main())
