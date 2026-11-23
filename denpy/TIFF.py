#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created 11 2026

Stacks of TIFF files to zarr or DEN, described by a JSON specification (docs/TIFFSTACK_SPEC.md).

	from denpy import TIFF
	plan = TIFF.validate("scan.json")         # or a dict; everything is checked before writing
	print(plan.describe())
	report = TIFF.execute(plan, threads=16)   # or TIFF.convert_spec("scan.json")

The specification has three parts, and this module one layer for each:
	tiffstacks     TiffStack       which files form a stack (reading)
	virtualstacks  LambdaStack,    stacks computed from other stacks: a lambda applied frame by
	               ReducedStack    frame, or a reduction of a stack to one frame (computing)
	arrays/array   ZarrWriter,     what is written where (writing)
	+ output       DenWriter
Tiff and virtual stacks offer the same interface (shape, dtype, frame(k), is_null(k)), and the
writers take anything with it.

Frames are computed and written in worker threads: a thread reads the TIFF frames of frame k,
applies the lambdas, converts and writes (zarr compresses in the same thread); zip stores are
written in parallel through ZARZS.ZipStore.

The command line interface is the script convertTiffStack.

@author: Vojtěch Kulvait
@license: GNU General Public License v3.0
"""
import ast
import collections
import copy
import importlib
import itertools
import json
import os
import shutil
import string
import tempfile
import threading
import time
from concurrent.futures import ThreadPoolExecutor

import numpy as np

SPEC = "denpy.tiffstack"
VERSION = 1
REDUCTIONS = ("mean", "median", "min", "max", "sum")
DEFAULT_MIN_THREADS, DEFAULT_MAX_THREADS = 16, 32


class SpecError(ValueError):
	"""The specification is invalid or does not match the files (found before anything is written)."""


class FrameError(RuntimeError):
	"""A frame does not match the specification while reading or converting it."""


def default_threads():
	"""A quarter of the CPUs available to the process, at least 16 and at most 32 (I/O bound work)."""
	try:
		cpus = len(os.sched_getaffinity(0))
	except AttributeError:
		cpus = os.cpu_count() or 1
	return min(DEFAULT_MAX_THREADS, max(DEFAULT_MIN_THREADS, cpus // 4))


def default_max_memory():
	"""A quarter of the RAM: the limit for reductions that hold a whole stack (median)."""
	try:
		return os.sysconf("SC_PAGE_SIZE") * os.sysconf("SC_PHYS_PAGES") // 4
	except (ValueError, OSError, AttributeError):
		return 8 << 30


# ---------------------------------------------------------------- helpers

def _keys(d, where, allowed, required=()):
	"""Check that d is an object with the required keys and no keys other than allowed ('_' comments)."""
	if not isinstance(d, dict):
		raise SpecError(f"{where}: expected an object, got {type(d).__name__}")
	unknown = [k for k in d if k not in allowed and not str(k).startswith("_")]
	if unknown:
		raise SpecError(f"{where}: unknown key(s) {', '.join(map(repr, unknown))}; allowed: {', '.join(allowed)}")
	missing = [k for k in required if k not in d]
	if missing:
		raise SpecError(f"{where}: missing key(s) {', '.join(map(repr, missing))}")


def _dtype(value, where):
	try:
		dtype = np.dtype(value)
	except TypeError as e:
		raise SpecError(f"{where}: invalid dtype {value!r}") from e
	if dtype.kind not in "biuf":
		raise SpecError(f"{where}: dtype {dtype} is not a number type")
	return dtype


def _int_list(value, where, length=None, minimum=0):
	if not isinstance(value, (list, tuple)) or not all(isinstance(v, int) and not isinstance(v, bool) for v in value):
		raise SpecError(f"{where}: expected a list of integers, got {value!r}")
	if length is not None and len(value) != length:
		raise SpecError(f"{where}: expected {length} integers, got {value!r}")
	if any(v < minimum for v in value):
		raise SpecError(f"{where}: values must be >= {minimum}, got {value!r}")
	return [int(v) for v in value]


def convert(array, dtype, clip=False, where="frame"):
	"""Convert to dtype: float to integer rounds to nearest; values outside the range of an integer
	dtype are clipped with clip=True, otherwise FrameError; NaN to integer is a FrameError."""
	array = np.asarray(array)
	if array.dtype == dtype:
		return array
	if dtype.kind in "iu" and array.size:
		if array.dtype.kind == "f":
			if np.isnan(array).any():
				raise FrameError(f"{where}: NaN can not be converted to {dtype}")
			array = np.rint(array)
		if array.dtype.kind in "biuf":
			info = np.iinfo(dtype)
			low, high = array.min(), array.max()
			if low < info.min or high > info.max:
				if not clip:
					raise FrameError(f"{where}: values in [{low}, {high}] outside the range of {dtype} "
									 f"[{info.min}, {info.max}]; set \"clip\": true to clip them")
				array = np.clip(array, info.min, info.max)
	return array.astype(dtype)


def compile_lambda(text, nargs, where):
	"""A function from a lambda expression string (checked with ast first: one lambda with nargs
	parameters; np is available) or a callable (python API)."""
	if callable(text):
		return text
	if not isinstance(text, str):
		raise SpecError(f"{where}: lambda must be a string, got {type(text).__name__}")
	try:
		tree = ast.parse(text.strip(), mode="eval")
	except SyntaxError as e:
		raise SpecError(f"{where}: invalid lambda {text!r}: {e.msg}") from e
	if not isinstance(tree.body, ast.Lambda):
		raise SpecError(f"{where}: expected a single lambda expression, got {text!r}")
	a = tree.body.args
	if a.vararg or a.kwarg or a.kwonlyargs or a.defaults:
		raise SpecError(f"{where}: the lambda must have plain parameters only: {text!r}")
	count = len(a.posonlyargs) + len(a.args)
	if count != nargs:
		raise SpecError(f"{where}: the lambda has {count} parameter(s), expected {nargs} "
						f"(the frame and {nargs - 1} from args): {text!r}")
	return eval(compile(tree, f"<{where}>", "eval"), {"np": np})


def import_function(name, where):
	"""The function module:name (name may be dotted, e.g. module:Class.method) or a callable."""
	if callable(name):
		return name
	if not isinstance(name, str) or name.count(":") != 1:
		raise SpecError(f"{where}: function must be 'module:name', got {name!r}")
	module, attribute = name.split(":")
	try:
		value = importlib.import_module(module)
		for part in attribute.split("."):
			value = getattr(value, part)
	except (ImportError, AttributeError) as e:
		raise SpecError(f"{where}: can not import {name!r}: {e}") from e
	if not callable(value):
		raise SpecError(f"{where}: {name!r} is not callable")
	return value


def _imap(executor, fn, items, window):
	"""executor.map with at most window tasks in flight; results in order."""
	pending = collections.deque()
	try:
		for item in items:
			pending.append(executor.submit(fn, item))
			if len(pending) >= window:
				yield pending.popleft().result()
		while pending:
			yield pending.popleft().result()
	finally:
		for future in pending:
			future.cancel()


class Stats:
	"""Counters of one phase (thread safe)."""

	def __init__(self):
		self._lock = threading.Lock()
		self.input_bytes = 0
		self.tiff_frames = 0

	def add_input(self, nbytes):
		with self._lock:
			self.input_bytes += nbytes
			self.tiff_frames += 1


# ---------------------------------------------------------------- reading

def read_tiff(path, page=0):
	"""A page of a TIFF file as a numpy array."""
	import tifffile
	with tifffile.TiffFile(path) as tif:
		try:
			return tif.pages[page].asarray()
		except IndexError:
			raise FrameError(f"{path} has no page {page}") from None


def tiff_header(path, page=0):
	"""(shape, dtype) of a page of a TIFF file without reading its data."""
	import tifffile
	with tifffile.TiffFile(path) as tif:
		try:
			p = tif.pages[page]
		except IndexError:
			raise FrameError(f"{path} has no page {page}") from None
		return tuple(p.shape), np.dtype(p.dtype)


def _pattern_files(entry, base_dir, where):
	_keys(entry, where, ("pattern", "range"), ("pattern", "range"))
	pattern = entry["pattern"]
	if not isinstance(pattern, str):
		raise SpecError(f"{where}: pattern must be a string")
	try:
		fields = [f for f in string.Formatter().parse(pattern) if f[1] is not None]
	except ValueError as e:
		raise SpecError(f"{where}: invalid pattern {pattern!r}: {e}") from e
	if len(fields) != 1 or fields[0][1] not in ("", "0"):
		raise SpecError(f"{where}: the pattern must have exactly one field {{}} or {{:...}} for the index: {pattern!r}")
	r = entry["range"]
	if not isinstance(r, list) or len(r) not in (2, 3) or not all(isinstance(v, int) and not isinstance(v, bool) for v in r):
		raise SpecError(f"{where}: range must be [start, stop] or [start, stop, step] (integers), got {r!r}")
	if len(r) == 3 and r[2] == 0:
		raise SpecError(f"{where}: range step must not be 0")
	try:
		return [(os.path.join(base_dir, pattern.format(i)), 0) for i in range(*r)]
	except (ValueError, TypeError) as e:
		raise SpecError(f"{where}: pattern {pattern!r} can not format an integer: {e}") from e


def expand_frames(entries, base_dir, where):
	"""Resolve the frame entries of a tiff stack: a list of (path, page) or None, and whether a
	pattern entry was used."""
	if not isinstance(entries, list):
		raise SpecError(f"{where}: frames must be a list")
	frames, has_pattern = [], False
	for i, entry in enumerate(entries):
		w = f"{where}[{i}]"
		if entry is None:
			frames.append(None)
		elif isinstance(entry, str):
			frames.append((os.path.join(base_dir, entry), 0))
		elif isinstance(entry, dict) and "pattern" in entry:
			frames.extend(_pattern_files(entry, base_dir, w))
			has_pattern = True
		elif isinstance(entry, dict):
			_keys(entry, w, ("file", "page"), ("file", ))
			page = entry.get("page", 0)
			if not isinstance(entry["file"], str) or not isinstance(page, int) or isinstance(page, bool) or page < 0:
				raise SpecError(f"{w}: expected {{\"file\": str, \"page\": int >= 0}}, got {entry!r}")
			frames.append((os.path.join(base_dir, entry["file"]), page))
		else:
			raise SpecError(f"{w}: a frame entry is a file name, null, {{\"file\": ..., \"page\": ...}} "
							f"or {{\"pattern\": ..., \"range\": ...}}, got {entry!r}")
	return frames, has_pattern


# ---------------------------------------------------------------- stacks

class Stack:
	"""Interface of tiff and virtual stacks: shape (n, Y, X) or (Y, X) for a single frame, dtype,
	frame(k) (k ignored for a single frame), is_null(k), fill (value of null frames or None)."""
	name = None
	shape = None
	dtype = None
	fill = None
	kind = None

	@property
	def nframes(self):
		"""Number of frames, None for a single frame (2D)."""
		return self.shape[0] if len(self.shape) == 3 else None

	@property
	def frame_shape(self):
		return tuple(self.shape[-2:])

	def is_null(self, k):
		return False

	def null_frames(self):
		return [k for k in range(self.nframes or 0) if self.is_null(k)]

	def dependencies(self):
		return []

	def prepare(self, executor, window, log):
		"""Compute what has to be in memory before frames are written (reductions)."""

	def provenance(self):
		raise NotImplementedError


class TiffStack(Stack):
	kind = "tiffstack"

	def __init__(self, name, spec, base_dir, check_headers=False):
		where = f"tiffstacks.{name}"
		_keys(spec, where, ("frames", "framecount", "frame_shape", "dtype", "fill_value"),
			  ("frames", "frame_shape", "dtype"))
		self.name, self.spec = name, spec
		self.frames, has_pattern = expand_frames(spec["frames"], base_dir, f"{where}.frames")
		if "framecount" in spec:
			count = spec["framecount"]
			if not isinstance(count, int) or isinstance(count, bool) or count != len(self.frames):
				raise SpecError(f"{where}: framecount is {count!r} but the frames resolve to {len(self.frames)} frames")
		elif has_pattern:
			raise SpecError(f"{where}: framecount is required with a pattern entry (it resolves to {len(self.frames)} frames)")
		frame_shape = _int_list(spec["frame_shape"], f"{where}.frame_shape", 2, minimum=1)
		self.dtype = _dtype(spec["dtype"], f"{where}.dtype")
		self.shape = (len(self.frames), *frame_shape)
		self._null = [k for k, f in enumerate(self.frames) if f is None]
		if self._null:
			if "fill_value" not in spec:
				raise SpecError(f"{where}: fill_value is required, the frames contain null (frames {self._null[:10]})")
			self.fill = _fill(spec["fill_value"], self.dtype, False, f"{where}.fill_value")
		elif "fill_value" in spec:
			self.fill = _fill(spec["fill_value"], self.dtype, False, f"{where}.fill_value")
		missing = [p for p in {f[0] for f in self.frames if f is not None} if not os.path.isfile(p)]
		if missing:
			raise SpecError(f"{where}: {len(missing)} listed file(s) do not exist, e.g. {sorted(missing)[:5]}")
		if check_headers:
			for f in self.frames:
				if f is not None:
					try:
						shape, dtype = tiff_header(*f)
					except FrameError as e:
						raise SpecError(f"{where}: {e}") from e
					if shape != tuple(frame_shape) or dtype != self.dtype:
						raise SpecError(f"{where}: {f[0]} page {f[1]} is {dtype} {list(shape)}, "
										f"the spec says {self.dtype} {frame_shape}")

	def is_null(self, k):
		return self.frames[k] is None

	def null_frames(self):
		return list(self._null)

	def frame(self, k, stats=None):
		entry = self.frames[k]
		if entry is None:
			return np.full(self.frame_shape, self.fill, dtype=self.dtype)
		path, page = entry
		try:
			data = read_tiff(path, page)
		except FrameError:
			raise
		except Exception as e:
			raise FrameError(f"tiffstacks.{self.name}[{k}]: can not read {path}: {e}") from e
		if data.shape != self.frame_shape or data.dtype != self.dtype:
			raise FrameError(f"tiffstacks.{self.name}[{k}]: {path} page {page} is {data.dtype} {list(data.shape)}, "
							 f"the spec says {self.dtype} {list(self.frame_shape)}")
		if stats is not None:
			stats.add_input(os.path.getsize(path) if page == 0 else data.nbytes)
		return data

	def provenance(self):
		frames = [None if f is None else (f[0] if f[1] == 0 else {"file": f[0], "page": f[1]}) for f in self.frames]
		spec = {k: v for k, v in self.spec.items() if k != "frames" and not k.startswith("_")}
		return {"kind": self.kind, **spec, "frames": frames}


def _fill(value, dtype, clip, where):
	if not isinstance(value, (int, float)) or isinstance(value, bool):
		raise SpecError(f"{where}: expected a number, got {value!r}")
	try:
		return convert(np.array(value), dtype, clip, where)[()]
	except FrameError as e:
		raise SpecError(str(e)) from None


VIRTUAL_KEYS = ("stack", "lambda", "function", "args", "reduce", "dtype", "clip", "fill_value")


class LambdaStack(Stack):
	"""Frame k: the lambda (identity without one) applied to frame k of the source and the args,
	converted to dtype."""
	kind = "lambda"

	def __init__(self, name, spec, source, args, where):
		self.name, self.spec, self.source, self.args, self.where = name, spec, source, args, where
		self.clip = spec.get("clip", False)
		if not isinstance(self.clip, bool):
			raise SpecError(f"{where}.clip: expected true or false")
		if "lambda" in spec and "function" in spec:
			raise SpecError(f"{where}: lambda and function exclude each other")
		if "lambda" in spec:
			self.function = compile_lambda(spec["lambda"], 1 + len(args), f"{where}.lambda")
		elif "function" in spec:
			self.function = import_function(spec["function"], f"{where}.function")
		else:
			if args:
				raise SpecError(f"{where}: args without lambda or function")
			self.function = None
		n = source.nframes
		for i, a in enumerate(args):
			if isinstance(a, Stack) and a.nframes is not None and a.nframes != n:
				raise SpecError(f"{where}.args[{i}]: {a.name} has {a.nframes} frames, "
								f"{source.name} {'has ' + str(n) if n is not None else 'is a single frame'}; "
								"an argument is a single frame (2D) or has as many frames as the stack")
		# the result frame: one run on frames of zeros
		if self.function is None:
			result_shape, result_dtype = source.frame_shape, source.dtype
		else:
			probe = [np.zeros(a.frame_shape, a.dtype) if isinstance(a, Stack) else a for a in args]
			try:
				with np.errstate(all="ignore"):
					result = np.asarray(self.function(np.zeros(source.frame_shape, source.dtype), *probe))
			except Exception as e:
				raise SpecError(f"{where}: the {'lambda' if 'lambda' in spec else 'function'} fails on frames of zeros "
								f"of the declared shapes and dtypes: {type(e).__name__}: {e}") from e
			if result.ndim != 2:
				raise SpecError(f"{where}: the result must be a 2D frame, got shape {list(result.shape)}")
			result_shape, result_dtype = result.shape, result.dtype
		self.dtype = _dtype(spec["dtype"], f"{where}.dtype") if "dtype" in spec else _dtype(result_dtype, where)
		self.shape = (n, *result_shape) if n is not None else tuple(result_shape)
		self._framewise = [a for a in [source, *args] if isinstance(a, Stack) and a.nframes is not None]
		self._null = sorted({k for a in self._framewise for k in a.null_frames()})
		fill = spec.get("fill_value")
		if fill is None:
			fill = next((a.fill for a in self._framewise if a.fill is not None), None)
			fill = None if fill is None else fill.item()
		if fill is not None:
			self.fill = _fill(fill, self.dtype, self.clip, f"{where}.fill_value")
		elif self._null:
			raise SpecError(f"{where}: frames {self._null[:10]} are null in a source without fill_value; set fill_value")
		self._single = None
		self._lock = threading.Lock()

	def dependencies(self):
		return [self.source, *[a for a in self.args if isinstance(a, Stack)]]

	def is_null(self, k):
		return any(a.is_null(k) for a in self._framewise)

	def null_frames(self):
		return list(self._null)

	def _compute(self, k, stats):
		x = self.source.frame(k, stats)
		if self.function is None:
			result = x
		else:
			values = [a.frame(k, stats) if isinstance(a, Stack) else a for a in self.args]
			result = np.asarray(self.function(x, *values))
			if result.shape != self.frame_shape:
				raise FrameError(f"{self.where}[{k}]: the result has shape {list(result.shape)}, "
								 f"expected {list(self.frame_shape)} (from the frames of zeros)")
		return convert(result, self.dtype, self.clip, f"{self.where}[{k}]")

	def frame(self, k=None, stats=None):
		if self.nframes is None:  # a single frame: computed once
			with self._lock:
				if self._single is None:
					self._single = self._compute(None, stats)
				return self._single
		if self.is_null(k):
			return np.full(self.frame_shape, self.fill, dtype=self.dtype)
		return self._compute(k, stats)

	def provenance(self):
		spec = {k: (v if not callable(v) else repr(v)) for k, v in self.spec.items() if not k.startswith("_")}
		return {"kind": "virtualstack", **spec}


class ReducedStack(Stack):
	"""One frame from a stack: mean, median, min, max or sum over its frames (null frames left out)."""
	kind = "reduce"

	def __init__(self, name, spec, source, where, max_memory):
		self.name, self.spec, self.source, self.where = name, spec, source, where
		self.op = spec["reduce"]
		if self.op not in REDUCTIONS:
			raise SpecError(f"{where}.reduce: expected one of {', '.join(REDUCTIONS)}, got {self.op!r}")
		for key in ("lambda", "function", "args"):
			if key in spec:
				raise SpecError(f"{where}: reduce and {key} exclude each other (reduce a virtual stack instead)")
		if source.nframes is None:
			raise SpecError(f"{where}: {source.name} is a single frame, nothing to reduce")
		self.frames = [k for k in range(source.nframes) if not source.is_null(k)]
		if not self.frames:
			raise SpecError(f"{where}: all frames of {source.name} are null")
		self.clip = spec.get("clip", False)
		if not isinstance(self.clip, bool):
			raise SpecError(f"{where}.clip: expected true or false")
		if "dtype" in spec:
			self.dtype = _dtype(spec["dtype"], f"{where}.dtype")
		elif self.op in ("min", "max"):
			self.dtype = source.dtype
		elif self.op == "sum":
			self.dtype = np.dtype(np.float64)
		else:
			self.dtype = np.dtype(np.float64 if source.dtype == np.float64 else np.float32)
		self.shape = source.frame_shape
		self.memory = len(self.frames) * int(np.prod(self.shape)) * source.dtype.itemsize if self.op == "median" else 0
		if self.memory > max_memory:
			raise SpecError(f"{where}: the median of {len(self.frames)} frames needs {self.memory / 2**30:.1f} GiB, "
							f"above the limit {max_memory / 2**30:.1f} GiB (max_memory)")
		self.value = None
		self.stats = None
		self.seconds = None

	def dependencies(self):
		return [self.source]

	def prepare(self, executor, window, log):
		if self.value is not None:
			return
		start, stats = time.time(), Stats()
		frames = _imap(executor, lambda k: self.source.frame(k, stats), self.frames, window)
		if self.op == "median":
			stack = np.empty((len(self.frames), *self.shape), dtype=self.source.dtype)
			for i, frame in enumerate(frames):
				stack[i] = frame
			result = np.median(stack, axis=0)
			del stack
		elif self.op in ("mean", "sum"):
			result = np.zeros(self.shape, dtype=np.float64)
			for frame in frames:
				result += frame
			if self.op == "mean":
				result /= len(self.frames)
		else:
			function = np.minimum if self.op == "min" else np.maximum
			result = None
			for frame in frames:
				result = frame.copy() if result is None else function(result, frame, out=result)
		self.value = convert(result, self.dtype, self.clip, self.where)
		self.stats, self.seconds = stats, time.time() - start
		log(f"TIME reduce {self.name} ({self.op} of {len(self.frames)} frames of {self.source.name}): {self.seconds:.2f}s")

	def frame(self, k=None, stats=None):
		if self.value is None:
			raise RuntimeError(f"{self.where} is used before it was computed")
		return self.value

	def provenance(self):
		return {"kind": "virtualstack", **{k: v for k, v in self.spec.items() if not k.startswith("_")}}


def make_virtual(name, spec, where, resolve, max_memory):
	"""A virtual stack from its spec; resolve(name) gives a source."""
	_keys(spec, where, VIRTUAL_KEYS, ("stack", ))
	source = resolve(spec["stack"], f"{where}.stack")
	if "reduce" in spec:
		return ReducedStack(name, spec, source, where, max_memory)
	args = spec.get("args", [])
	if not isinstance(args, list):
		raise SpecError(f"{where}.args: expected a list")
	values = []
	for i, a in enumerate(args):
		if isinstance(a, str):
			values.append(resolve(a, f"{where}.args[{i}]"))
		elif isinstance(a, (int, float)) and not isinstance(a, bool):
			values.append(a)
		else:
			raise SpecError(f"{where}.args[{i}]: expected a source name or a number, got {a!r}")
	return LambdaStack(name, spec, source, values, where)


# ---------------------------------------------------------------- plan

ARRAY_KEYS = VIRTUAL_KEYS + ("grid", "attributes", "overwrite", "zarr")
ZARR_KEYS = ("codec", "chunks", "shards", "dimension_names")
OUTPUT_KEYS = ("format", "path", "store", "mode", "zarr_format", "attributes", "atomic", "provenance")


class OutputArray:
	"""An array to write: a stack, the grid of its frames and the storage settings."""

	def __init__(self, name, spec, stack, output):
		where = "array" if name is None else f"arrays.{name}"
		self.name, self.spec, self.stack, self.where = name, spec, stack, where
		n = stack.nframes
		if n is None:
			if "grid" in spec:
				raise SpecError(f"{where}.grid: {stack.name or name} is a single frame (2D), no grid")
			self.grid = []
		else:
			self.grid = _int_list(spec.get("grid", [n]), f"{where}.grid", minimum=0 if n == 0 else 1)
			if not self.grid or int(np.prod(self.grid)) != n:
				raise SpecError(f"{where}.grid: {self.grid} has {int(np.prod(self.grid)) if self.grid else 0} frames, "
								f"the stack has {n}")
		self.shape = (*self.grid, *stack.frame_shape)
		self.dtype = stack.dtype
		self.attributes = spec.get("attributes", {})
		if not isinstance(self.attributes, dict):
			raise SpecError(f"{where}.attributes: expected an object")
		self.overwrite = spec.get("overwrite", False)
		if not isinstance(self.overwrite, bool):
			raise SpecError(f"{where}.overwrite: expected true or false")
		if output["format"] == "zarr":
			self._zarr(spec.get("zarr", {}), output, f"{where}.zarr")
		else:  # zarr settings and attributes are ignored, so one spec can be written to either format
			from denpy import DEN
			try:
				DEN.npDtypeToDenDataType(self.dtype)
			except KeyError:
				raise SpecError(f"{where}: DEN can not store dtype {self.dtype}") from None

	def _zarr(self, settings, output, where):
		_keys(settings, where, ZARR_KEYS)
		ndim = len(self.shape)
		self.codec = settings.get("codec", {"name": "none"})
		if not isinstance(self.codec, dict) or not isinstance(self.codec.get("name"), str):
			raise SpecError(f"{where}.codec: expected {{\"name\": ..., \"clevel\": ..., <codec arguments>}}, got {self.codec!r}")
		chunks = settings.get("chunks", [1] * len(self.grid) + [-1, -1])
		self.chunks = [self.shape[i] if c == -1 else c for i, c in enumerate(_int_list(chunks, f"{where}.chunks", ndim, minimum=-1))]
		if any(c == 0 for c in self.chunks):
			raise SpecError(f"{where}.chunks: 0 is not a chunk size")
		self.shards = None
		v2 = output.get("zarr_format", 3) == 2
		if "shards" in settings:
			if v2:
				raise SpecError(f"{where}.shards: only zarr_format 3")
			shards = _int_list(settings["shards"], f"{where}.shards", ndim, minimum=-1)
			self.shards = [self.shape[i] if s == -1 else s for i, s in enumerate(shards)]
			if any(s == 0 or s % c for s, c in zip(self.shards, self.chunks)):
				raise SpecError(f"{where}.shards: {self.shards} must be multiples of chunks {self.chunks}")
		self.dimension_names = settings.get("dimension_names")
		if self.dimension_names is not None:
			if v2:
				raise SpecError(f"{where}.dimension_names: only zarr_format 3")
			if not isinstance(self.dimension_names, list) or len(self.dimension_names) != ndim:
				raise SpecError(f"{where}.dimension_names: expected {ndim} names")
		# frames are written in blocks along the grid axes that never share a chunk / shard
		self.block = (self.shards or self.chunks)[:len(self.grid)]

	def blocks(self):
		"""Index blocks along the grid axes (tuples of slices), each written at once; a block never
		shares a chunk (shard) with another, so blocks can be written in parallel."""
		if not self.grid:
			return [()]
		block = getattr(self, "block", [1] * len(self.grid))
		starts = itertools.product(*[range(0, g, b) for g, b in zip(self.grid, block)])
		return [tuple(slice(s, min(s + b, g)) for s, b, g in zip(start, block, self.grid)) for start in starts]

	def frames_of(self, block):
		"""Frame indices k (C order of the grid) of a block, in the C order of the block."""
		if not self.grid:
			return [None]
		return [int(np.ravel_multi_index(index, self.grid))
				for index in itertools.product(*[range(s.start, s.stop) for s in block])]


class Plan:
	"""A validated specification: sources, arrays to write, output settings."""

	def __init__(self, spec, base_dir, sources, arrays, output, single):
		self.spec, self.base_dir, self.sources, self.arrays, self.output, self.single = \
			spec, base_dir, sources, arrays, output, single

	def needed(self):
		"""Sources needed by the arrays, dependencies first."""
		order, seen = [], set()

		def visit(stack):
			if id(stack) in seen:
				return
			seen.add(id(stack))
			for d in stack.dependencies():
				visit(d)
			order.append(stack)
		for a in self.arrays:
			visit(a.stack)
		return order

	def describe(self):
		lines = [f"output: {self.output['format']} {self.output['path']}"
				 + (f" ({self.output['store']} store, zarr v{self.output.get('zarr_format', 3)})" if self.output["format"] == "zarr" else "")
				 + f", mode {self.output.get('mode', 'w-')}"]
		for s in self.needed():
			nulls = s.null_frames()
			lines.append(f"  {s.kind:9s} {s.name or '':20s} {str(list(s.shape)):22s} {s.dtype}"
						 + (f"  null frames {nulls[:10]}{' ...' if len(nulls) > 10 else ''}" if nulls else ""))
		for a in self.arrays:
			extra = f"  chunks {a.chunks} codec {a.codec}" if self.output["format"] == "zarr" else ""
			lines.append(f"  write     {a.name or '(root)':20s} {str(list(a.shape)):22s} {a.dtype}{extra}")
		return "\n".join(lines)


def load_spec(spec, base_dir=None):
	"""(spec dict, base_dir) from a JSON file name or a dict (copied)."""
	if isinstance(spec, (str, os.PathLike)):
		path = os.fspath(spec)
		with open(path) as f:
			try:
				data = json.load(f)
			except json.JSONDecodeError as e:
				raise SpecError(f"{path}: invalid JSON: {e}") from e
		spec_dir = os.path.dirname(os.path.abspath(path))
	elif isinstance(spec, dict):
		data, spec_dir = copy.copy(spec), os.getcwd()
	else:
		raise SpecError(f"expected a file name or a dict, got {type(spec).__name__}")
	if base_dir is None:
		base_dir = os.path.join(spec_dir, data["base_dir"]) if isinstance(data.get("base_dir"), str) else spec_dir
	return data, base_dir


def _merge(defaults, own):
	merged = dict(defaults)
	for key, value in own.items():
		if isinstance(value, dict) and isinstance(merged.get(key), dict):
			merged[key] = {**merged[key], **value}
		else:
			merged[key] = value
	return merged


def _check_output(output, single):
	_keys(output, "output", OUTPUT_KEYS, ("format", "path"))
	fmt = output["format"]
	if fmt not in ("zarr", "den"):
		raise SpecError(f"output.format: expected \"zarr\" or \"den\", got {fmt!r}")
	if not isinstance(output["path"], str) or not output["path"]:
		raise SpecError("output.path: expected a path")
	mode = output.setdefault("mode", "w-")
	if mode not in ("w-", "w", "a"):
		raise SpecError(f"output.mode: expected \"w-\", \"w\" or \"a\", got {mode!r}")
	if mode == "a" and single:
		raise SpecError("output.mode \"a\" adds arrays to a group / directory: use arrays, not array")
	for key in ("atomic", "provenance"):
		if not isinstance(output.setdefault(key, key == "atomic"), bool):
			raise SpecError(f"output.{key}: expected true or false")
	if fmt == "zarr":
		if output.get("store") not in ("zip", "directory"):
			raise SpecError(f"output.store: expected \"zip\" or \"directory\", got {output.get('store')!r}")
		if output.setdefault("zarr_format", 3) not in (2, 3):
			raise SpecError("output.zarr_format: expected 2 or 3")
		if not isinstance(output.setdefault("attributes", {}), dict):
			raise SpecError("output.attributes: expected an object")
		if single and output["attributes"]:
			raise SpecError("output.attributes: the root is the array, use the attributes of the array")
	else:
		for key in ("store", "zarr_format", "attributes"):
			if key in output:
				raise SpecError(f"output.{key}: only for zarr output")
		if output["provenance"]:
			raise SpecError("output.provenance: only for zarr output")
	path = output["path"]
	if fmt == "den" and not single:
		# a directory of DEN files is a container, never replaced as a whole: the files are
		# checked one by one (_check_existing)
		if os.path.lexists(path) and not os.path.isdir(path):
			raise SpecError(f"output.path: {path} exists and is no directory")
		if mode == "a" and not os.path.isdir(path):
			raise SpecError(f"output.path: {path} does not exist (mode \"a\")")
		return
	if mode == "w-" and os.path.lexists(path):
		raise SpecError(f"output.path: {path} exists (mode \"w-\"); use mode \"w\" to replace it")
	if mode == "a" and not os.path.exists(path):
		raise SpecError(f"output.path: {path} does not exist (mode \"a\")")


def _valid_name(name, where):
	parts = name.split("/") if isinstance(name, str) else [""]
	if any(p in ("", ".", "..") for p in parts):
		raise SpecError(f"{where}: invalid name {name!r}")


def validate(spec, base_dir=None, check_headers=False, max_memory=None):
	"""A Plan from a specification (file name or dict); SpecError if anything does not fit."""
	spec, base_dir = load_spec(spec, base_dir)
	_keys(spec, "spec", ("spec", "version", "base_dir", "tiffstacks", "virtualstacks", "array", "arrays",
						 "defaults", "output"), ("spec", "version", "tiffstacks", "output"))
	if spec["spec"] != SPEC:
		raise SpecError(f"spec: expected {SPEC!r}, got {spec['spec']!r}")
	if spec["version"] != VERSION:
		raise SpecError(f"version: only version {VERSION} is supported, got {spec['version']!r}")
	if ("array" in spec) == ("arrays" in spec):
		raise SpecError("exactly one of array and arrays is required")
	single = "array" in spec
	if single and "defaults" in spec:
		raise SpecError("defaults: only with arrays")
	max_memory = default_max_memory() if max_memory is None else max_memory
	output = dict(spec["output"]) if isinstance(spec["output"], dict) else spec["output"]
	_check_output(output, single)

	tiff_specs = spec["tiffstacks"]
	virtual_specs = spec.get("virtualstacks", {})
	for key, value in (("tiffstacks", tiff_specs), ("virtualstacks", virtual_specs)):
		if not isinstance(value, dict):
			raise SpecError(f"{key}: expected an object name -> stack")
	both = set(tiff_specs) & set(virtual_specs)
	if both:
		raise SpecError(f"names used in both tiffstacks and virtualstacks: {sorted(both)}")
	sources, building = {}, []

	def resolve(name, where):
		if not isinstance(name, str):
			raise SpecError(f"{where}: expected the name of a source, got {name!r}")
		if name in sources:
			return sources[name]
		if name in building:
			raise SpecError(f"{where}: cycle {' -> '.join(building[building.index(name):] + [name])}")
		if name in tiff_specs:
			sources[name] = TiffStack(name, tiff_specs[name], base_dir, check_headers)
		elif name in virtual_specs:
			building.append(name)
			sources[name] = make_virtual(name, virtual_specs[name], f"virtualstacks.{name}", resolve, max_memory)
			building.pop()
		else:
			raise SpecError(f"{where}: no tiff stack or virtual stack named {name!r} "
							f"(arrays can not be referred to; define it as a virtual stack)")
		return sources[name]

	for name in list(tiff_specs) + list(virtual_specs):
		resolve(name, name)

	if single:
		array_specs = {None: spec["array"]}
	else:
		if not isinstance(spec["arrays"], dict) or not spec["arrays"]:
			raise SpecError("arrays: expected an object name -> array with at least one array")
		defaults = spec.get("defaults", {})
		_keys(defaults, "defaults", ARRAY_KEYS)
		array_specs = {name: _merge(defaults, a) if isinstance(a, dict) else a for name, a in spec["arrays"].items()}
	arrays = []
	for name, a in array_specs.items():
		where = "array" if name is None else f"arrays.{name}"
		if name is not None:
			_valid_name(name, where)
		_keys(a, where, ARRAY_KEYS, ("stack", ))
		computed = {k: v for k, v in a.items() if k in VIRTUAL_KEYS}
		if set(computed) == {"stack"}:
			stack = resolve(a["stack"], f"{where}.stack")
		else:
			stack = make_virtual(name or "array", computed, where, resolve, max_memory)
		arrays.append(OutputArray(name, a, stack, output))
	if not single:
		_check_existing(arrays, output)
	return Plan(spec, base_dir, sources, arrays, output, single)


def _check_existing(arrays, output):
	"""Arrays that exist in the output: mode "a" needs overwrite for them; for a directory of DEN
	files also mode "w-" refuses them ("w" replaces them)."""
	if output["format"] == "den":
		exists = [a for a in arrays if os.path.lexists(os.path.join(output["path"], a.name + ".den"))]
		if output["mode"] == "w-" and exists:
			raise SpecError(f"{[a.name + '.den' for a in exists]} exist in {output['path']} (mode \"w-\"); "
							"use mode \"w\" to replace them")
		if output["mode"] == "w":
			return
	elif output["mode"] != "a":
		return
	else:
		import zarr
		store = _open_store(output, output["path"], "r")
		try:
			group = zarr.open_group(store, mode="r")
			exists = [a for a in arrays if a.name in group]
		except Exception as e:
			raise SpecError(f"output.path: can not open {output['path']} as a zarr group: {e}") from e
		finally:
			store.close()
	clash = [a.name for a in exists if not a.overwrite]
	if clash:
		raise SpecError(f"arrays {clash} exist in {output['path']}; set \"overwrite\": true to replace them")


# ---------------------------------------------------------------- writing

def _open_store(output, path, mode):
	import zarr
	if output["store"] == "zip":
		from denpy import ZARZS
		return ZARZS.ZipStore(path, mode=mode)
	return zarr.storage.LocalStore(path, read_only=(mode == "r"))


def _tree_bytes(path):
	if os.path.isfile(path):
		return os.path.getsize(path)
	return sum(os.path.getsize(os.path.join(d, f)) for d, _, files in os.walk(path) for f in files)


class _Target:
	"""Where the output goes: with atomic writing a temporary name next to path, renamed at the end."""

	def __init__(self, output, is_dir):
		self.output, self.path, self.is_dir = output, output["path"], is_dir
		self.atomic = output["atomic"] and output["mode"] != "a"
		if not self.atomic:
			self.work = self.path
			if output["mode"] == "w" and os.path.lexists(self.path):
				_remove(self.path)
			if is_dir:
				os.makedirs(self.path, exist_ok=True)
		else:
			parent = os.path.dirname(os.path.abspath(self.path))
			prefix = os.path.basename(os.path.normpath(self.path)) + ".tmp-"
			if is_dir:
				self.work = tempfile.mkdtemp(prefix=prefix, dir=parent)
			else:
				fd, self.work = tempfile.mkstemp(prefix=prefix, dir=parent)
				os.close(fd)

	def commit(self):
		if not self.atomic:
			return
		if os.path.lexists(self.path):
			if self.output["mode"] == "w-":
				raise FileExistsError(f"{self.path} appeared while writing; the result is in {self.work}")
			_remove(self.path)
		os.replace(self.work, self.path)

	def abort(self):
		if self.atomic and os.path.lexists(self.work):
			_remove(self.work)


def _remove(path):
	if os.path.isdir(path) and not os.path.islink(path):
		shutil.rmtree(path)
	else:
		os.remove(path)


class ZarrWriter:

	def __init__(self, plan):
		import zarr
		self.plan, self.output = plan, plan.output
		self.target = _Target(self.output, is_dir=self.output["store"] == "directory")
		mode = self.output["mode"]
		self.store = _open_store(self.output, self.target.work, "a" if mode == "a" else "w")
		self.zarr_format = self.output["zarr_format"]
		self.group = None
		if not plan.single:
			self.group = zarr.open_group(self.store, mode="a" if mode == "a" else "w", zarr_format=self.zarr_format)
			if self.output["attributes"]:
				self.group.attrs.update(self.output["attributes"])

	def create(self, a):
		import zarr
		from denpy import ZAR
		kwargs = dict(shape=a.shape, dtype=a.dtype, chunks=a.chunks, attributes=a.attributes or None)
		if a.stack.fill is not None:
			kwargs["fill_value"] = a.stack.fill
		if a.shards:
			kwargs["shards"] = a.shards
		if a.dimension_names:
			kwargs["dimension_names"] = a.dimension_names
		codec = dict(a.codec)
		name = codec.pop("name")
		if name == "none":
			kwargs["compressors"] = None
		else:
			clevel = codec.pop("clevel", 5)
			kwargs.update(ZAR.codec_arguments(ZAR.get_compressor(name, clevel, zarrv2=self.zarr_format == 2, dtype=a.dtype, **codec)))
		if a.name is None:
			return zarr.create_array(self.store, zarr_format=self.zarr_format, overwrite=True, **kwargs)
		return self.group.create_array(a.name, overwrite=a.overwrite or self.output["mode"] != "a", **kwargs)

	def write(self, array, block, data):
		array[block] = data

	def set_attributes(self, array, attributes):
		array.attrs.update(attributes)

	def stored_bytes(self, a):
		"""Bytes of the array in the output (after finish)."""
		path = self.output["path"]
		prefix = "" if a.name is None else a.name + "/"
		if self.output["store"] == "zip":
			import zipfile
			with zipfile.ZipFile(path) as z:
				return sum(i.compress_size for i in z.infolist() if i.filename.startswith(prefix))
		return _tree_bytes(os.path.join(path, a.name) if a.name else path)

	def finish(self):
		self.store.close()
		self.target.commit()

	def abort(self):
		try:
			self.store.close()
		finally:
			self.target.abort()


class DenWriter:

	"""One DEN file (array), or DEN files <name>.den in a directory (arrays). The directory is a
	container: it is created if needed and never removed; each file is written under a temporary
	name next to it and renamed at the end, other files of the directory are not touched."""

	def __init__(self, plan):
		self.plan, self.output = plan, plan.output
		self.target = _Target(self.output, is_dir=False) if plan.single else None
		self.fds = []

	def file(self, a):
		return self.target.work if a.name is None else os.path.join(self.output["path"], a.name + ".den")

	def create(self, a):
		from denpy import DEN
		path = self.file(a)
		os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
		work = path
		if a.name is not None:  # each file of a directory through a temporary name
			fd, work = tempfile.mkstemp(prefix=os.path.basename(path) + ".tmp-", dir=os.path.dirname(os.path.abspath(path)))
			os.close(fd)
		DEN.writeEmptyDEN(work, list(reversed(a.shape)), elementtype=a.dtype, force=True)
		fd = os.open(work, os.O_RDWR)
		self.fds.append((fd, work, path))
		frame_bytes = int(np.prod(a.stack.frame_shape)) * a.dtype.itemsize
		return {"fd": fd, "frame_bytes": frame_bytes, "dtype": DEN.get_dtype_with_byte_order(a.dtype), "array": a}

	def write(self, handle, block, data):
		data = np.ascontiguousarray(data, dtype=handle["dtype"])
		a = handle["array"]
		frames = a.frames_of(block)
		for i, k in enumerate(frames):  # frames of a block are not consecutive in general
			frame = data.reshape(-1, *a.stack.frame_shape)[i]
			offset = 4096 + (k or 0) * handle["frame_bytes"]
			view, written = memoryview(frame).cast("B"), 0
			while written < view.nbytes:
				written += os.pwrite(handle["fd"], view[written:], offset + written)

	def set_attributes(self, handle, attributes):
		pass

	def stored_bytes(self, a):
		return os.path.getsize(os.path.join(self.output["path"], a.name + ".den") if a.name else self.output["path"])

	def finish(self):
		for fd, work, path in self.fds:
			os.close(fd)
		for fd, work, path in self.fds:
			if work != path:
				if self.output["mode"] == "w-" and os.path.lexists(path):
					raise FileExistsError(f"{path} appeared while writing; the result is in {work}")
				os.replace(work, path)
		self.fds = []
		if self.target is not None:
			self.target.commit()

	def abort(self):
		for fd, work, path in self.fds:
			try:
				os.close(fd)
			except OSError:
				pass
			if work != path and os.path.exists(work):
				os.remove(work)
		self.fds = []
		if self.target is not None:
			self.target.abort()


def execute(plan, threads=None, verbose=False, log=print):
	"""Write the arrays of a Plan (or of a specification, validated first); returns the report."""
	if not isinstance(plan, Plan):
		plan = validate(plan)
	threads = default_threads() if threads is None or threads <= 0 else threads
	window = 2 * threads
	start = time.time()
	report = {"output": plan.output["path"], "format": plan.output["format"], "threads": threads,
			  "virtualstacks": {}, "arrays": {}}
	writer = (ZarrWriter if plan.output["format"] == "zarr" else DenWriter)(plan)
	try:
		with ThreadPoolExecutor(threads) as executor:
			for source in plan.needed():
				source.prepare(executor, window, log)
				if isinstance(source, ReducedStack):
					report["virtualstacks"][source.name] = {"reduce": source.op, "frames": len(source.frames),
															"input_bytes": source.stats.input_bytes,
															"seconds": round(source.seconds, 3)}
			handles = []
			for a in plan.arrays:
				t0, stats = time.time(), Stats()
				handle = writer.create(a)
				handles.append(handle)

				def task(block, a=a, handle=handle, stats=stats):
					frames = [a.stack.frame(k, stats) for k in a.frames_of(block)]
					data = frames[0] if not a.grid else np.stack(frames).reshape(
						*[s.stop - s.start for s in block], *a.stack.frame_shape)
					writer.write(handle, block, data)
					if verbose:
						log(f"  {a.name or 'array'} {[(s.start, s.stop) for s in block]}")

				for _ in _imap(executor, task, a.blocks(), window):
					pass
				seconds = time.time() - t0
				nulls = a.stack.null_frames()
				report["arrays"][a.name or ""] = {"shape": list(a.shape), "dtype": str(a.dtype),
												  "frames": a.stack.nframes or 1, "null_frames": nulls,
												  "tiff_frames_read": stats.tiff_frames,
												  "input_bytes": stats.input_bytes, "seconds": round(seconds, 3)}
				log(f"TIME write {a.name or 'array'} {list(a.shape)} {a.dtype}: {seconds:.2f}s")
			if plan.output["format"] == "zarr" and plan.output["provenance"]:
				for a, handle in zip(plan.arrays, handles):
					writer.set_attributes(handle, {"denpy_tiffstack": _provenance(plan, a, report["arrays"][a.name or ""])})
		writer.finish()
	except BaseException:
		writer.abort()
		raise
	for a in plan.arrays:
		entry = report["arrays"][a.name or ""]
		entry["stored_bytes"] = writer.stored_bytes(a)
		raw = int(np.prod(a.shape)) * a.dtype.itemsize
		entry["ratio_to_raw"] = round(entry["stored_bytes"] / raw, 4) if raw else None
	report["seconds"] = round(time.time() - start, 3)
	log(f"TIME total: {report['seconds']:.2f}s")
	return report


def _provenance(plan, a, array_report):
	"""The array spec, the sources it is computed from (tiff stacks with expanded frames) and its report."""
	named = any(a.stack is s for s in plan.sources.values())
	sources, todo = {}, [a.stack] if named else list(a.stack.dependencies())
	while todo:
		s = todo.pop()
		if s.name not in sources:
			sources[s.name] = s.provenance()
			todo.extend(s.dependencies())
	spec = {k: (v if not callable(v) else repr(v)) for k, v in a.spec.items() if not k.startswith("_")}
	return {"spec": SPEC, "version": VERSION, "array": spec, "sources": sources, "report": array_report}


def convert_spec(spec, threads=None, check_headers=False, max_memory=None, verbose=False, log=print):
	"""validate and execute in one call; returns the report."""
	return execute(validate(spec, check_headers=check_headers, max_memory=max_memory), threads, verbose, log)


def spec_from_files(files, output, frame_shape=None, dtype=None, base_dir=None, **array):
	"""A spec of one tiff stack written as one array (the root of a zarr store or one DEN file);
	frame_shape and dtype from the first file if not given."""
	files = [os.fspath(f) for f in files]
	if not files:
		raise SpecError("no files")
	if frame_shape is None or dtype is None:
		shape, first_dtype = tiff_header(os.path.join(base_dir or os.getcwd(), files[0]))
		frame_shape = list(shape) if frame_shape is None else frame_shape
		dtype = str(first_dtype) if dtype is None else dtype
	spec = {"spec": SPEC, "version": VERSION,
			"tiffstacks": {"tiffs": {"frames": files, "framecount": len(files),
									 "frame_shape": list(frame_shape), "dtype": str(np.dtype(dtype))}},
			"array": {"stack": "tiffs", **array},
			"output": output}
	if base_dir is not None:
		spec["base_dir"] = os.fspath(base_dir)
	return spec


EXTENSION_KEYS = ("tiffstacks", "virtualstacks", "arrays", "defaults")


def extend_spec(spec, extension):
	"""A copy of spec extended by extension (a dict or a JSON file name), e.g. a spec generated from
	metadata extended by processing:

		tiffstacks     new tiff stacks (an existing name is an error)
		virtualstacks  added; an existing name is replaced
		arrays         added; an existing name is replaced; null removes the array
		defaults       merged into the defaults (dicts such as zarr key by key)

	The result is validated as any spec (validate)."""
	if isinstance(extension, (str, os.PathLike)):
		with open(extension) as f:
			try:
				extension = json.load(f)
			except json.JSONDecodeError as e:
				raise SpecError(f"{os.fspath(extension)}: invalid JSON: {e}") from e
	_keys(extension, "extension", EXTENSION_KEYS)
	result = copy.deepcopy(spec)
	for key in ("tiffstacks", "virtualstacks", "arrays", "defaults"):
		if key in extension and not isinstance(extension[key], dict):
			raise SpecError(f"extension.{key}: expected an object")
	clash = sorted(set(extension.get("tiffstacks", {})) & set(result.get("tiffstacks", {})))
	if clash:
		raise SpecError(f"extension.tiffstacks: {clash} exist in the spec")
	result.setdefault("tiffstacks", {}).update(copy.deepcopy(extension.get("tiffstacks", {})))
	if extension.get("virtualstacks"):
		result.setdefault("virtualstacks", {}).update(copy.deepcopy(extension["virtualstacks"]))
	if "arrays" in extension:
		if "array" in result:
			raise SpecError("extension.arrays: the spec has a single array (array), not arrays")
		arrays = result.setdefault("arrays", {})
		for name, value in extension["arrays"].items():
			if value is None:
				if name not in arrays:
					raise SpecError(f"extension.arrays.{name}: null removes an array, there is no such array")
				del arrays[name]
			else:
				arrays[name] = copy.deepcopy(value)
	if extension.get("defaults"):
		if "array" in result:
			raise SpecError("extension.defaults: the spec has a single array (array), no defaults")
		result["defaults"] = _merge(result.get("defaults", {}), copy.deepcopy(extension["defaults"]))
	return result
