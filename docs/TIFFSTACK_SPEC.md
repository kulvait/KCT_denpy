# TIFF stack specification (version 1)

Implemented by `denpy/TIFF.py` (tests: `denpy/test_TIFF.py`) and the script `convertTiffStack`.

A JSON document describing how stacks of TIFF files become arrays in zarr or DEN. It has three
separate parts, and `denpy.TIFF` has one layer for each:

| part | layer | what it says |
|---|---|---|
| `tiffstacks` | reading | which files form a stack, their frame shape and dtype, which frames are known to be missing |
| `virtualstacks` | computing | stacks computed from other stacks (a lambda applied frame by frame, or a reduction to one frame); not written |
| `arrays` / `array` | writing | what is written: a stack as it is or computed the same way as a virtual stack, with the storage settings |
| `output` | writing | where and how the arrays are stored: zarr (zip or directory store) or DEN |

The shape and dtype of every array follow from the spec alone, before any TIFF is read, so a spec
can be validated completely up front; the files are then checked against it while reading.

The specification says *what* to build; execution options (threads, verbosity, report file)
belong to the CLI `convertTiffStack` / the function call, not to the spec. Metadata parsers (P05
HDF5, nanoCT log, ...) produce such a spec and add their metadata to the output separately (zarr
`mode: "a"`), so the converter knows nothing about beamlines.

All ranges are half-open, `[start, stop)`, as in python: `[0, 20]` is 0..19.
Arrays are described in numpy order: a frame is `[Y, X]`, a 3D array `[frames, Y, X]`.

## Examples

`...` only shortens this document; a real spec lists every frame or uses a pattern entry.

One stack written as it is; the zarr root is the array itself:

```json
{
  "spec": "denpy.tiffstack", "version": 1,
  "tiffstacks": {
    "proj": {"frames": [{"pattern": "img{}.tif", "range": [2, 501]}], "framecount": 499,
             "frame_shape": [2048, 2048], "dtype": "uint16"}
  },
  "array": {"stack": "proj"},
  "output": {"format": "zarr", "store": "zip", "path": "proj.zar"}
}
```

A P05 / nanoCT style scan with the usual names `dar`, `fla`, `img`, written as it is (a root group
with three arrays):

```json
{
  "spec": "denpy.tiffstack", "version": 1,
  "base_dir": "/gpfs/.../raw/scan_042",
  "tiffstacks": {
    "dar": {"frames": ["dar_00000.tif", "dar_00001.tif", ..., "dar_00019.tif"],
            "frame_shape": [3840, 5120], "dtype": "uint16"},
    "fla": {"frames": [{"pattern": "ref_{:05d}.tif", "range": [0, 20]},
                       {"pattern": "ref_{:05d}.tif", "range": [1520, 1540]}],
            "framecount": 40, "frame_shape": [3840, 5120], "dtype": "uint16"},
    "img": {"frames": ["img_00000.tif", "img_00001.tif", null, "img_00003.tif", ..., "img_01499.tif"],
            "framecount": 1500, "frame_shape": [3840, 5120], "dtype": "uint16", "fill_value": 0}
  },
  "arrays": {
    "dar": {"stack": "dar"},
    "fla": {"stack": "fla"},
    "img": {"stack": "img", "zarr": {"dimension_names": ["angle", "y", "x"]}}
  },
  "defaults": {"zarr": {"codec": {"name": "htj2k", "clevel": 0}}},
  "output": {"format": "zarr", "store": "zip", "path": "scan_042.zar", "mode": "w"}
}
```

`null`: the metadata lists frame 2 but there is no file for it; the spec says so explicitly and
`fill_value` is written there, the report lists the filled frames.

The same scan dark corrected, into DEN (`scan_042_den/dar.den`, `dar_median.den`, `fla.den`,
`img.den`):

```json
{
  "spec": "denpy.tiffstack", "version": 1,
  "tiffstacks": {"dar": {...}, "fla": {...}, "img": {...}},
  "virtualstacks": {
    "dar_median": {"stack": "dar", "reduce": "median", "dtype": "uint16"}
  },
  "arrays": {
    "dar": {"stack": "dar"},
    "dar_median": {"stack": "dar_median"},
    "fla": {"stack": "fla", "dtype": "uint16", "clip": true,
            "lambda": "lambda x, dark: x.astype('int32') - dark.astype('int32')",
            "args": ["dar_median"]},
    "img": {"stack": "img", "dtype": "uint16", "clip": true,
            "lambda": "lambda x, dark: x.astype('int32') - dark.astype('int32')",
            "args": ["dar_median"]}
  },
  "output": {"format": "den", "path": "scan_042_den"}
}
```

`dar_median` is computed once; leaving out its entry in `arrays` uses it without writing it.

## Top level

| key | required | meaning |
|---|---|---|
| `spec` | yes | `"denpy.tiffstack"` |
| `version` | yes | `1` |
| `base_dir` | no | relative TIFF paths resolve against it; default: the directory of the spec file (cwd for a spec passed as a dict) |
| `tiffstacks` | yes | object name → tiff stack, see Tiff stacks |
| `virtualstacks` | no | object name → virtual stack, see Virtual stacks |
| `array` | one of | a single array, see Arrays: the zarr root is this array; DEN: `output.path` is the file |
| `arrays` | one of | object name → array: the zarr root is a group with these arrays, names are paths in it (`"raw/img"` creates group `raw`); DEN: `output.path` is a directory with `<name>.den` |
| `defaults` | no | only with `arrays`: any array key; an array's own keys override it (dicts such as `zarr` merge key by key) |
| `output` | yes | see Output |

Exactly one of `array` and `arrays`.

Names: tiff stacks and virtual stacks are the *sources* and share one namespace (a name is used
once in both). `stack` and `args` refer to sources only. Arrays are outputs: their names are paths
in the output, may repeat a source name (`"img": {"stack": "img"}`) and can not be referred to;
to write a computed stack and use it elsewhere, define it as a virtual stack and write it with
`{"stack": name}` (it is computed once). Sources that no array needs are neither read nor computed.

Unknown keys are an error (typos must not pass silently); keys starting with `_` are ignored
(comments).

## Tiff stacks

| key | default | meaning |
|---|---|---|
| `frames` | required | list of frame entries, see below |
| `framecount` | | number of frames the entries resolve to; **required** when `frames` contains a pattern entry, optional otherwise; when given it must match, checked at validation |
| `frame_shape` | required | `[Y, X]` of every frame |
| `dtype` | required | dtype of every frame (`"uint16"`, `"float32"`, ...) |
| `fill_value` | | value of the `null` frames; **required** when `frames` contains `null` |

`frames` is a list of entries, in stack order. An entry is:

| entry | frames | meaning |
|---|---|---|
| `"a.tif"` | 1 | page 0 of a TIFF file; relative to `base_dir` |
| `{"file": "a.tif", "page": 3}` | 1 | one page of a multi-page TIFF |
| `null` | 1 | a frame that should exist but has no file (e.g. listed in the metadata, absent on disk); `fill_value` |
| `{"pattern": "ref_{:05d}.tif", "range": [start, stop]}` | `stop - start` | the files `pattern.format(i)` for `i` in `range(start, stop)`; `"range": [start, stop, step]` as python's `range` |

`pattern` is python `str.format` with exactly one integer field (`{}`, `{:05d}`, ...), nothing
else is interpreted. Because a pattern hides the count, a stack with a pattern entry must state
`framecount`; validation expands every entry and compares, so a wrong range is caught before
anything is written. The expanded list is what the report and `provenance` record.

Everything a stack lists must be there: validation checks that all files exist (and with
`--check-headers` that their TIFF shape and dtype match `frame_shape` and `dtype`); while reading,
a frame that can not be read or does not match is an error. The only frames not read from a file
are the `null` entries, which the spec author put there knowingly.

## Virtual stacks

A virtual stack says how to get frame k from other sources. It is one of

| form | shape | meaning |
|---|---|---|
| `{"stack": "img", "lambda": "...", "args": [...]}` | `[n, Y', X']` | the lambda applied to every frame |
| `{"stack": "dark", "reduce": "median"}` | `[Y, X]` | one frame from the whole stack: `"mean"`, `"median"`, `"min"`, `"max"`, `"sum"` over the frames |

| key | default | meaning |
|---|---|---|
| `stack` | required | the source frame k comes from (`n` = its number of frames) |
| `lambda` | | a python lambda expression, see Lambdas |
| `function` | | instead of `lambda`: `"module:name"`, an importable function with the same signature |
| `args` | `[]` | further arguments of the lambda after the frame: a source name or a number |
| `reduce` | | `"mean"`, `"median"`, `"min"`, `"max"` or `"sum"` (not together with `lambda`) |
| `dtype` | of the result | dtype of the frames; the result is converted to it, see Conversion. Reductions without `dtype`: `mean` and `median` float32 (float64 for a float64 stack), `min` and `max` the dtype of the stack, `sum` float64 |
| `clip` | `false` | `true`: values outside the range of an integer `dtype` are set to its minimum / maximum; `false`: they are an error |
| `fill_value` | of the sources | value of null frames (see Conversion); default: the `fill_value` of `stack`, else of the first stack in `args` that has one |

References form a graph that must not have cycles.

### Lambdas

```
"lambda": "lambda x, dark: x.astype('int32') - dark.astype('int32')", "args": ["dar_median"]
```

The lambda gets frame k of `stack` and then one value per entry of `args`:

- a 2D source (a reduction, shape `[Y, X]`): that frame, the same for every k
- a source with `n` frames: its frame k (frame-by-frame operations between stacks, e.g. two cameras);
  any other number of frames is an error
- a number: the number

It returns a 2D array. `np` (numpy) is available; the string must be a single `lambda` expression
with plain parameters, 1 + `len(args)` of them, checked with `ast` at validation before anything is
evaluated. Examples:

```
"lambda x: x[100:3940, :]"                                     crop
"lambda x: x.reshape(1920, 2, 2560, 2).mean(axis=(1, 3))"       2x2 binning
"lambda x, dark, flat: (x - dark) / np.maximum(flat - dark, 1)" flat field (args: two reductions)
"lambda x, k: k * x"                                             args: ["gain"] or a number, e.g. [0.5]
```

Validation runs each lambda once on frames of zeros of the declared shapes and dtypes to learn the
shape and dtype of its result; while writing, every frame must give that same shape. A spec with
lambdas is code: run only specs you trust (as for any script).

`function: "module:name"` (also `module:Class.method`) is the same with an importable function, for
anything longer than a line, e.g. a hot pixel filter; in the python API `lambda` and `function` may
also be a callable defined in the executing code. The function is called from several threads at
once.

### Conversion

The result of a lambda or reduction is converted to `dtype`: float to integer rounds to nearest
(`np.rint`); values outside the integer range are clipped with `clip: true`, otherwise an error
naming the frame. A frame of a virtual stack is null when it is null in `stack` or in a stack of
`args` with `n` frames: it is not given to the lambda, the result frame is `fill_value` converted
to `dtype`. Reductions leave null frames out (a stack of only `null` frames can not be reduced).

## Arrays

An array is what is written. It is a source as it is, `{"stack": "img"}`, or it is computed like a
virtual stack (the same keys `stack`, `lambda`/`function`, `args`, `reduce`, `dtype`, `clip`),
which is the same as a virtual stack written with `{"stack": ...}`. Further keys:

| key | default | meaning |
|---|---|---|
| `grid` | `[n]` | shape of the frame grid in numpy order, product = `n`; frames fill it in C order (last index fastest). `[30, 50]` → array shape `(30, 50, Y, X)` |
| `attributes` | `{}` | zarr attributes of the array |
| `overwrite` | `false` | with output `mode: "a"`: replace an existing array (zarr) or file (DEN) of this name |
| `zarr` | | `codec`, `chunks`, `shards`, `dimension_names`, see Output |

DEN output ignores `attributes` and `zarr`, so one spec can be written to either format.

## Output

| key | formats | default | meaning |
|---|---|---|---|
| `format` | | required | `"zarr"` or `"den"` |
| `path` | both | required | zarr: the store (any name, the suffix is not interpreted; `.zar` recommended for zip stores); DEN: the file (`array`) or directory (`arrays`) |
| `store` | zarr | required | `"zip"` (one file, `ZARZS.ZipStore`) or `"directory"` |
| `mode` | both | `"w-"` | `"w-"` fail if the output exists, `"w"` replace it, `"a"` add arrays to an existing zarr group / DEN directory (only with `arrays`; an existing array of the same name is an error unless the array sets `"overwrite": true`) |
| `zarr_format` | zarr | `3` | `2` or `3` |
| `attributes` | zarr | `{}` | attributes of the root group (with `arrays`; with `array` use the array's `attributes`) |
| `atomic` | both | `true` | write to `path.tmp-XXXX` and rename at the end (only for `mode` `"w"`/`"w-"`) |
| `provenance` | zarr | `false` | attribute `denpy_tiffstack` of each array: its spec, the sources it is computed from (tiff stacks with the expanded frame lists, `null` included) and its report |

Per array `zarr`:

| key | default | meaning |
|---|---|---|
| `codec` | `{"name": "none"}` | `{"name": ..., "clevel": ..., <codec kwargs>}`, i.e. the arguments of `ZAR.get_compressor`; names as in `--compression` (`none`, `zstd`, `lz4`, `blosc-*`, `gzip`, `htj2k`, `jpeg2k`, `jpegxl`, `avif`, ...) |
| `chunks` | `[1, ..., 1, -1, -1]` | one frame per chunk; `-1` = full extent of that axis. Frames are written in blocks of the chunk (shard) size along the frame axes, so larger chunks along the frames work in parallel too |
| `shards` | none | zarr v3 sharding, same notation (multiple of `chunks`) |
| `dimension_names` | none | zarr v3 |

A DEN directory (`arrays` with format `den`) is a container: it may exist and hold other files
(e.g. a working directory), it is created if needed and never removed. Each `<name>.den` is written
under a temporary name next to it and renamed when all arrays are written; `"w-"` refuses existing
files of the arrays, `"w"` replaces just them, `"a"` needs `overwrite` for them. A failed conversion
leaves the existing files as they were.

For DEN, `/` in an array name becomes a subdirectory (`raw/img` → `<path>/raw/img.den`). DEN
writes an array with dimspec = numpy shape reversed: `(X, Y, grid[-1], ..., grid[0])`, the
existing DEN extended header; DEN dtypes are uint8, uint16, int16, uint32, int32, uint64, int64,
float32, float64 (others are an error at validation).

## Execution

1. Validate: schema, names, references (no cycles), expand frame lists, check `framecount`, files
   exist, lambdas (syntax, parameter count, one run on zero frames), shape and dtype of every
   source and array, frame counts of `args`, DEN dtypes, memory of reductions (see below).
2. Compute the reductions needed, in dependency order (kept in memory).
3. Write the arrays: worker threads compute frame k (read the TIFF frames it needs, apply the
   lambdas, convert) and write it (compression in the same thread); with ZARZS the zip is written
   in parallel too. Frame k of an intermediate virtual stack is computed when needed, not stored.

Reductions: `mean`, `sum`, `min`, `max` accumulate frame by frame; `median` needs the whole stack in
memory (20 darks of 5120×3840 float32: 1.6 GB) and validation rejects it above a limit
(`--max-memory`, default a quarter of the RAM).

The layers are separate in the code: tiff stacks and virtual stacks both offer `shape`, `dtype` and
`frame(k)`, and the writers take anything with these, so a later version can read stacks from DEN
or zarr files as well.

## CLI

```
convertTiffStack spec.json [-j THREADS] [--report report.json] [--dry-run] [--check-headers]
                 [--max-memory GiB] [--verbose]
convertTiffStack --files img*.tif --zarr out.zar [--store zip|directory] [--zarr-format 2|3]
                 [--compression htj2k --clevel 0] [--dtype float32 [--clip]] [--sort natural]
                 [--mode w-|w] [--write-spec spec.json] [--dry-run]
convertTiffStack --files img*.tif --den out.den ...
```

The second form builds a spec with one tiff stack and one `array` from the file list (in the order
given, so the shell glob order; `--sort natural` sorts `img2` before `img10`; `frame_shape` and
`dtype` from the first file) and runs it; `--write-spec` saves that spec as a starting point for
editing. `--dry-run` validates and prints the sources and arrays with shapes, dtypes and null
frames. Threads default to a quarter of the CPUs, at least 16, at most 32. Exit code 1 with an
`ERROR:` line for an invalid spec or a frame that does not match.

## Extensions

`TIFF.extend_spec(spec, extension)` (CLI: `convertTiffStack spec.json --extend EXT.json`) returns the
spec extended by a JSON object with the keys

| key | effect |
|---|---|
| `tiffstacks` | new tiff stacks (an existing name is an error) |
| `virtualstacks` | added; an existing name is replaced |
| `arrays` | added; an existing name is replaced; `null` removes the array |
| `defaults` | merged into the defaults (dicts such as `zarr` key by key) |

This is how a spec generated from metadata gets processing added, e.g. a dark median and dark
corrected stacks.

## From scan datasets

denpy has no measurement specific code here: the specs of tomography scans are built in lautools,
`TOMOGRAPHY.scanTiffStackConversionSpecs(scanDataset, output, ...)` from the DataFrame of
`PETRA.scanDataset` or `NANOCT.scanDataset`, then run with `denpy.TIFF`. The lautools command
`tiffScanDataToZarrForMicrotomography` is the worked example.

## Report

Execution returns (and with `--report out.json` writes):

- per array: `shape`, `dtype`, `frames`, `null_frames`, `tiff_frames_read`, `input_bytes` (TIFF file
  sizes read), `stored_bytes`, `ratio_to_raw` (stored / uncompressed array), `seconds`
- per reduction: `reduce`, `frames`, `input_bytes`, `seconds`
- `threads`, total `seconds`

The CLI prints `TIME` lines per reduction and array and a `SIZE` line per array.

## Not in version 1

- pyramids (run `ZARPYRAMID` on the result)
- binning along the stack (frame k from several frames of the source)
- stacks from DEN or zarr files (the layers allow it)
