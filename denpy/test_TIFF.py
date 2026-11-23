#!/usr/bin/env python
#-*- coding: utf-8 -*-
"""
Created 2026

Test of TIFF: TIFF stack specifications to zarr and DEN, validation, virtual stacks, convertTiffStack.

@author: Vojtech Kulvait
"""
import contextlib
import io
import json
import os
import tempfile
import unittest
import zipfile

import numpy as np

try:
	import tifffile
	import zarr
except ImportError:  # optional dependencies (denpy[zarr])
	tifffile = zarr = None

from denpy import DEN

if tifffile is not None:
	from denpy import TIFF, ZAR
	from denpy.scripts import convertTiffStack

Y, X = 24, 32
SUBTRACT = "lambda x, dark: x.astype('int32') - dark.astype('int32')"


def quiet(message):
	pass


@unittest.skipIf(tifffile is None, "tifffile and zarr are needed")
class TestTiffStack(unittest.TestCase):

	def setUp(self):
		self.tmp = tempfile.TemporaryDirectory()
		self.dir = self.tmp.name
		rng = np.random.default_rng(1)
		self.data = {}
		for prefix, n, base in (("dar", 5, 100), ("fla", 6, 3000), ("img", 10, 1500)):
			frames = (base + rng.integers(0, 200, (n, Y, X))).astype(np.uint16)
			for i, frame in enumerate(frames):
				tifffile.imwrite(self.path(f"{prefix}_{i:05d}.tif"), frame)
			self.data[prefix] = frames
		os.remove(self.path("img_00004.tif"))  # listed in the metadata, absent on disk
		self.median = np.rint(np.median(self.data["dar"], axis=0)).astype(np.uint16)
		ZAR.register_codecs()

	def tearDown(self):
		self.tmp.cleanup()

	def path(self, name):
		return os.path.join(self.dir, name)

	def spec(self, output, **extra):
		spec = {
			"spec": "denpy.tiffstack", "version": 1, "base_dir": self.dir,
			"tiffstacks": {
				"dar": {"frames": [{"pattern": "dar_{:05d}.tif", "range": [0, 5]}], "framecount": 5,
						"frame_shape": [Y, X], "dtype": "uint16"},
				"fla": {"frames": [f"fla_{i:05d}.tif" for i in range(6)], "frame_shape": [Y, X], "dtype": "uint16"},
				"img": {"frames": [f"img_{i:05d}.tif" if i != 4 else None for i in range(10)], "framecount": 10,
						"frame_shape": [Y, X], "dtype": "uint16", "fill_value": 0}},
			"virtualstacks": {"dar_median": {"stack": "dar", "reduce": "median", "dtype": "uint16"}},
			"arrays": {
				"dar": {"stack": "dar"},
				"dar_median": {"stack": "dar_median"},
				"fla": {"stack": "fla", "dtype": "uint16", "clip": True, "lambda": SUBTRACT, "args": ["dar_median"]},
				"img": {"stack": "img", "dtype": "uint16", "clip": True, "lambda": SUBTRACT, "args": ["dar_median"]}},
			"output": output}
		spec.update(extra)
		return spec

	def corrected(self, name):
		result = np.clip(self.data[name].astype(np.int32) - self.median, 0, 65535).astype(np.uint16)
		if name == "img":
			result[4] = 0
		return result

	def run_spec(self, spec, **kwargs):
		return TIFF.convert_spec(spec, threads=4, log=quiet, **kwargs)

	# ---- writing

	def check_scan(self, read):
		np.testing.assert_array_equal(read("dar"), self.data["dar"])
		np.testing.assert_array_equal(read("dar_median"), self.median)
		np.testing.assert_array_equal(read("fla"), self.corrected("fla"))
		np.testing.assert_array_equal(read("img"), self.corrected("img"))

	def test_scan_to_zarr_zip(self):
		zar = self.path("scan.zar")
		report = self.run_spec(self.spec({"format": "zarr", "store": "zip", "path": zar, "provenance": True},
										 defaults={"zarr": {"codec": {"name": "htj2k", "clevel": 0}}}))
		self.assertEqual(report["arrays"]["img"]["null_frames"], [4])
		self.assertEqual(report["arrays"]["img"]["tiff_frames_read"], 9)
		with zipfile.ZipFile(zar) as z:
			self.assertIsNone(z.testzip())
		group = zarr.open_group(zarr.storage.ZipStore(zar, mode="r"), mode="r")
		self.check_scan(lambda name: group[name][:])
		provenance = group["img"].attrs["denpy_tiffstack"]
		self.assertEqual(sorted(provenance["sources"]), ["dar", "dar_median", "img"])
		self.assertIsNone(provenance["sources"]["img"]["frames"][4])
		self.assertEqual(provenance["report"]["null_frames"], [4])
		self.assertFalse([f for f in os.listdir(self.dir) if ".tmp-" in f])

	def test_scan_to_zarr_directory_v2(self):
		out = self.path("scan.zarr")
		self.run_spec(self.spec({"format": "zarr", "store": "directory", "path": out, "zarr_format": 2},
								defaults={"zarr": {"codec": {"name": "zstd", "clevel": 3}}}))
		group = zarr.open_group(out, mode="r")
		self.assertEqual(group.metadata.zarr_format, 2)
		self.check_scan(lambda name: group[name][:])

	def test_scan_to_den(self):
		self.run_spec(self.spec({"format": "den", "path": self.path("den")}))
		self.assertEqual(sorted(os.listdir(self.path("den"))), ["dar.den", "dar_median.den", "fla.den", "img.den"])
		self.check_scan(lambda name: DEN.getNumpyArray(self.path(f"den/{name}.den")))
		self.assertEqual(DEN.readHeader(self.path("den/img.den"))["type"], np.uint16)

	def test_den_directory_is_a_container(self):
		# DEN files into an existing directory (e.g. the working directory): the other files stay, the
		# directory is never replaced; "w-" refuses existing arrays, "w" replaces just them
		wd = self.path("wd")
		os.makedirs(os.path.join(wd, "log"))
		for name in ("notes.txt", "img.den", "log/run.out"):
			with open(os.path.join(wd, name), "w") as f:
				f.write("keep")
		spec = self.spec({"format": "den", "path": wd})
		with self.assertRaisesRegex(TIFF.SpecError, "img.den.*exist"):
			TIFF.validate(spec)
		spec["output"]["mode"] = "w"
		self.run_spec(spec)
		self.assertEqual(sorted(os.listdir(wd)), ["dar.den", "dar_median.den", "fla.den", "img.den", "log", "notes.txt"])
		with open(os.path.join(wd, "log/run.out")) as f:
			self.assertEqual(f.read(), "keep")
		np.testing.assert_array_equal(DEN.getNumpyArray(os.path.join(wd, "img.den")), self.corrected("img"))
		# a failure leaves the existing files as they were and no temporary files
		before = {f: os.path.getmtime(os.path.join(wd, f)) for f in os.listdir(wd)}
		spec["arrays"]["img"]["clip"] = False
		spec["arrays"]["img"]["lambda"] = "lambda x, dark: x.astype('int32') - dark.astype('int32') - 5000"
		with self.assertRaises(TIFF.FrameError):
			self.run_spec(spec)
		self.assertEqual({f: os.path.getmtime(os.path.join(wd, f)) for f in os.listdir(wd)}, before)

	def test_single_array_and_grid(self):
		# the root of the zarr store is the array; a 4D array from a grid of frames
		spec = self.spec({"format": "zarr", "store": "zip", "path": self.path("fla.zar")})
		del spec["arrays"]
		spec["array"] = {"stack": "fla", "grid": [2, 3], "dtype": "float32", "zarr": {"chunks": [1, 3, -1, -1]},
						 "attributes": {"kind": "flat"}}
		self.run_spec(spec)
		array = zarr.open_array(zarr.storage.ZipStore(self.path("fla.zar"), mode="r"), mode="r")
		self.assertEqual(array.shape, (2, 3, Y, X))
		self.assertEqual(array.chunks, (1, 3, Y, X))
		self.assertEqual(array.attrs["kind"], "flat")
		np.testing.assert_array_equal(array[:], self.data["fla"].astype(np.float32).reshape(2, 3, Y, X))
		spec["output"] = {"format": "den", "path": self.path("fla.den")}
		del spec["array"]["zarr"], spec["array"]["attributes"]
		self.run_spec(spec)
		self.assertEqual(DEN.readHeader(self.path("fla.den"))["dimspec"], (X, Y, 3, 2))
		np.testing.assert_array_equal(DEN.getNumpyArray(self.path("fla.den")), self.data["fla"].reshape(2, 3, Y, X))

	def test_lambdas_functions_reductions(self):
		spec = self.spec({"format": "den", "path": self.path("den")})
		spec["virtualstacks"].update({
			"dar_mean": {"stack": "dar", "reduce": "mean"},
			"fla_norm": {"stack": "fla", "lambda": "lambda x, dark: x - dark", "args": ["dar_mean"]}})
		spec["arrays"] = {
			"crop": {"stack": "img", "lambda": "lambda x: x[2:20, 4:28]"},
			"bin": {"stack": "fla", "lambda": f"lambda x: x.reshape({Y // 2}, 2, {X // 2}, 2).mean(axis=(1, 3))"},
			"flip": {"stack": "dar", "function": "numpy:fliplr"},
			"ffc": {"stack": "img", "lambda": "lambda x, dark, flat, s: s * (x - dark) / np.maximum(flat, 1)",
					"args": ["dar_mean", "fla_mean", 2.0], "dtype": "float32"},
			"pairs": {"stack": "dar", "lambda": "lambda x, y: x.astype('float32') - y", "args": ["fla5"]},
			"dar_min": {"stack": "dar", "reduce": "min"},
			"dar_max": {"stack": "dar", "reduce": "max"},
			"dar_sum": {"stack": "dar", "reduce": "sum"},
			"fla_mean": {"stack": "fla_norm", "reduce": "mean"}}
		spec["virtualstacks"]["fla_mean"] = {"stack": "fla_norm", "reduce": "mean"}
		spec["virtualstacks"]["fla5"] = {"stack": "fla", "lambda": "lambda x: x[:]"}
		spec["tiffstacks"]["fla"]["frames"] = spec["tiffstacks"]["fla"]["frames"][:5]
		self.run_spec(spec)
		read = lambda name: DEN.getNumpyArray(self.path(f"den/{name}.den"))
		dar, fla, img = (self.data[k].astype(np.float64) for k in ("dar", "fla", "img"))
		fla = fla[:5]
		dar_mean = dar.mean(axis=0).astype(np.float32)
		fla_mean = (fla.astype(np.float32) - dar_mean).mean(axis=0, dtype=np.float64).astype(np.float32)
		self.assertEqual(read("crop").shape, (10, 18, 24))
		np.testing.assert_array_equal(read("crop")[0], self.data["img"][0, 2:20, 4:28])
		np.testing.assert_allclose(read("bin")[1], fla[1].reshape(Y // 2, 2, X // 2, 2).mean(axis=(1, 3)))
		np.testing.assert_array_equal(read("flip"), self.data["dar"][:, :, ::-1])
		expected = 2.0 * (img[0] - dar_mean) / np.maximum(fla_mean, 1)
		np.testing.assert_allclose(read("ffc")[0], expected, rtol=1e-5)
		np.testing.assert_array_equal(read("ffc")[4], 0)
		np.testing.assert_array_equal(read("pairs"), dar.astype(np.float32) - fla.astype(np.float32))
		np.testing.assert_array_equal(read("dar_min"), self.data["dar"].min(axis=0))
		np.testing.assert_array_equal(read("dar_max"), self.data["dar"].max(axis=0))
		np.testing.assert_array_equal(read("dar_sum"), dar.sum(axis=0))
		self.assertEqual(read("dar_sum").dtype, np.float64)

	def test_pages_of_a_multipage_tiff(self):
		tifffile.imwrite(self.path("stack.tif"), self.data["fla"])
		spec = {"spec": "denpy.tiffstack", "version": 1, "base_dir": self.dir,
				"tiffstacks": {"s": {"frames": [{"file": "stack.tif", "page": p} for p in (5, 0, 3)],
									 "frame_shape": [Y, X], "dtype": "uint16"}},
				"array": {"stack": "s"}, "output": {"format": "den", "path": self.path("s.den")}}
		self.run_spec(spec)
		np.testing.assert_array_equal(DEN.getNumpyArray(self.path("s.den")), self.data["fla"][[5, 0, 3]])

	def test_append_to_zip(self):
		zar = self.path("scan.zar")
		spec = self.spec({"format": "zarr", "store": "zip", "path": zar})
		spec["arrays"] = {"dar": {"stack": "dar"}}
		self.run_spec(spec)
		spec["output"]["mode"] = "a"
		spec["arrays"] = {"img": {"stack": "img"}, "dar": {"stack": "dar_median"}}
		with self.assertRaisesRegex(TIFF.SpecError, "overwrite"):
			TIFF.validate(spec)
		spec["arrays"]["dar"]["overwrite"] = True
		self.run_spec(spec)
		group = zarr.open_group(zarr.storage.ZipStore(zar, mode="r"), mode="r")
		np.testing.assert_array_equal(group["dar"][:], self.median)
		self.assertEqual(group["img"].shape, (10, Y, X))
		with zipfile.ZipFile(zar) as z:
			self.assertIsNone(z.testzip())
			self.assertEqual(len(z.namelist()), len(set(z.namelist())))

	def test_failure_leaves_no_output(self):
		# uint16 minus the dark without clip: negative values are an error, nothing is written
		spec = self.spec({"format": "zarr", "store": "zip", "path": self.path("scan.zar")})
		spec["arrays"]["img"]["clip"] = False
		spec["arrays"]["img"]["lambda"] = "lambda x, dark: x.astype('int32') - dark.astype('int32') - 1600"
		with self.assertRaisesRegex(TIFF.FrameError, "clip"):
			self.run_spec(spec)
		self.assertEqual(sorted(f for f in os.listdir(self.dir) if not f.endswith(".tif")), [])

	# ---- validation

	def assertSpecError(self, spec, pattern, **kwargs):
		with self.assertRaisesRegex(TIFF.SpecError, pattern):
			TIFF.validate(spec, **kwargs)

	def test_validation(self):
		out = {"format": "den", "path": self.path("den")}
		cases = []

		def case(pattern, change):
			spec = json.loads(json.dumps(self.spec(dict(out))))
			change(spec)
			cases.append((spec, pattern))

		case("framecount is required", lambda s: s["tiffstacks"]["dar"].pop("framecount"))
		case("framecount is 6", lambda s: s["tiffstacks"]["dar"].update(framecount=6))
		case("framecount is 7 but the frames resolve to 6", lambda s: s["tiffstacks"]["fla"].update(framecount=7))
		case("do not exist", lambda s: s["tiffstacks"]["fla"]["frames"].append("fla_00099.tif"))
		case("fill_value is required", lambda s: s["tiffstacks"]["img"].pop("fill_value"))
		case("exactly one field", lambda s: s["tiffstacks"]["dar"]["frames"][0].update(pattern="dar_{}_{}.tif"))
		case("unknown key.*'frame_count'", lambda s: s["tiffstacks"]["dar"].update(frame_count=5))
		case("unknown key.*'lamda'", lambda s: s["arrays"]["img"].update(lamda="lambda x: x"))
		case("2 parameter\\(s\\), expected 1", lambda s: s["arrays"]["img"].pop("args"))
		case("single lambda expression", lambda s: s["arrays"]["img"].update(**{"lambda": "x + 1", "args": []}))
		case("fails on frames of zeros", lambda s: s["arrays"]["img"].update(**{"lambda": "lambda x, d: x @ d"}))
		case("2D frame", lambda s: s["arrays"]["img"].update(**{"lambda": "lambda x, d: x.ravel()"}))
		case("arrays can not be referred to", lambda s: s["arrays"]["img"].update(args=["fla"]) or s["tiffstacks"].pop("fla"))
		case("has 6 frames", lambda s: s["arrays"]["img"].update(args=["fla"]))
		case("cycle", lambda s: s["virtualstacks"].update(a={"stack": "b"}, b={"stack": "a", "lambda": "lambda x: x"}))
		case("both tiffstacks and virtualstacks", lambda s: s["virtualstacks"].update(dar={"stack": "fla", "reduce": "mean"}))
		case("DEN can not store dtype int8", lambda s: s["arrays"]["img"].update(dtype="int8"))
		case("grid", lambda s: s["arrays"]["img"].update(grid=[3, 3]))
		case("exactly one of array and arrays", lambda s: s.update(array={"stack": "dar"}))
		case("reduce and lambda", lambda s: s["virtualstacks"]["dar_median"].update(**{"lambda": "lambda x: x"}))
		for spec, pattern in cases:
			with self.subTest(pattern=pattern):
				self.assertSpecError(spec, pattern)
		self.assertSpecError(self.spec(dict(out)), "needs 0.0 GiB, above the limit", max_memory=100)
		os.makedirs(self.path("den"))
		TIFF.validate(self.spec(dict(out)))  # an existing directory is fine, it is a container
		open(self.path("den/img.den"), "w").close()
		self.assertSpecError(self.spec(dict(out)), "img.den.*exist")

	def test_check_headers(self):
		spec = self.spec({"format": "den", "path": self.path("den")})
		spec["tiffstacks"]["fla"]["dtype"] = "float32"
		TIFF.validate(spec)  # only the existence of the files
		self.assertSpecError(spec, "the spec says float32", check_headers=True)
		with self.assertRaisesRegex(TIFF.FrameError, "the spec says float32"):
			self.run_spec(spec)

	# ---- extensions

	def test_extend_spec(self):
		base = self.spec({"format": "den", "path": self.path("den")})
		del base["virtualstacks"], base["arrays"]["dar_median"], base["arrays"]["fla"]
		base["arrays"]["img"] = {"stack": "img"}
		base["defaults"] = {"zarr": {"codec": {"name": "zstd"}}, "dtype": "uint16"}
		extension = {"_comment": "dark correction",
					 "virtualstacks": {"dar_median": {"stack": "dar", "reduce": "median", "dtype": "uint16"}},
					 "arrays": {"dar": None, "dar_median": {"stack": "dar_median"},
								"img": {"stack": "img", "clip": True, "lambda": SUBTRACT, "args": ["dar_median"]}},
					 "defaults": {"zarr": {"chunks": [1, -1, -1]}}}
		spec = TIFF.extend_spec(base, extension)
		self.assertEqual(sorted(spec["arrays"]), ["dar_median", "img"])
		self.assertEqual(spec["defaults"], {"zarr": {"codec": {"name": "zstd"}, "chunks": [1, -1, -1]}, "dtype": "uint16"})
		self.assertIn("dar", base["arrays"])  # the base spec is not changed
		self.run_spec(spec)
		np.testing.assert_array_equal(DEN.getNumpyArray(self.path("den/img.den")), self.corrected("img"))
		self.assertEqual(sorted(os.listdir(self.path("den"))), ["dar_median.den", "img.den"])
		for bad, pattern in (({"output": {}}, "unknown key"), ({"arrays": {"x": None}}, "no such array"),
							 ({"tiffstacks": {"dar": {}}}, "exist in the spec"), ({"arrays": []}, "expected an object")):
			with self.subTest(pattern=pattern):
				with self.assertRaisesRegex(TIFF.SpecError, pattern):
					TIFF.extend_spec(base, bad)

	# ---- convertTiffStack

	def cli(self, *argv):
		out = io.StringIO()
		with contextlib.redirect_stdout(out), contextlib.redirect_stderr(out):
			code = convertTiffStack.main([str(a) for a in argv])
		return code, out.getvalue()

	def test_cli(self):
		files = [self.path(f"dar_{i:05d}.tif") for i in (3, 1, 4, 0, 2)]
		code, out = self.cli("--files", *files, "--sort", "natural", "--zarr", self.path("dar.zar"),
							 "--compression", "htj2k", "--clevel", "0", "--write-spec", self.path("dar.json"),
							 "--report", self.path("report.json"))
		self.assertEqual(code, 0, out)
		self.assertIn("SIZE array", out)
		array = zarr.open_array(zarr.storage.ZipStore(self.path("dar.zar"), mode="r"), mode="r")
		np.testing.assert_array_equal(array[:], self.data["dar"])
		with open(self.path("report.json")) as f:
			self.assertEqual(json.load(f)["arrays"][""]["shape"], [5, Y, X])
		# the saved spec runs again (to DEN)
		with open(self.path("dar.json")) as f:
			spec = json.load(f)
		spec["output"] = {"format": "den", "path": self.path("dar.den")}
		with open(self.path("dar.json"), "w") as f:
			json.dump(spec, f)
		code, out = self.cli(self.path("dar.json"), "--dry-run")
		self.assertEqual(code, 0, out)
		self.assertFalse(os.path.exists(self.path("dar.den")))
		self.assertEqual(self.cli(self.path("dar.json"))[0], 0)
		np.testing.assert_array_equal(DEN.getNumpyArray(self.path("dar.den")), self.data["dar"])
		code, out = self.cli(self.path("dar.json"))
		self.assertEqual(code, 1)
		self.assertIn("exists", out)
		# a spec file with an extension, relative paths stay relative to the spec file
		spec["output"] = {"format": "den", "path": self.path("ext_den")}
		spec.pop("base_dir")
		spec["tiffstacks"]["tiffs"]["frames"] = [os.path.basename(p) for p in spec["tiffstacks"]["tiffs"]["frames"]]
		spec["arrays"] = {"dar": {"stack": "tiffs"}}
		del spec["array"]
		with open(self.path("rel.json"), "w") as f:
			json.dump(spec, f)
		with open(self.path("ext.json"), "w") as f:
			json.dump({"virtualstacks": {"m": {"stack": "tiffs", "reduce": "max"}}, "arrays": {"dar_max": {"stack": "m"}}}, f)
		cwd = os.getcwd()
		try:
			os.chdir("/")
			code, out = self.cli(self.path("rel.json"), "--extend", self.path("ext.json"), "--write-spec", self.path("full.json"))
		finally:
			os.chdir(cwd)
		self.assertEqual(code, 0, out)
		np.testing.assert_array_equal(DEN.getNumpyArray(self.path("ext_den/dar_max.den")), self.data["dar"].max(axis=0))
		with open(self.path("full.json")) as f:
			self.assertIn("m", json.load(f)["virtualstacks"])


if __name__ == "__main__":
	unittest.main()
