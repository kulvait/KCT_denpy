# denpy

**denpy** is a Python package designed for manipulating the [DEN file format](https://kulvait.github.io/KCT_doc/den-format.html). The main module `DEN` contains functions for conversion between `DEN` and `numpy.array`. Other modules include functions and methods for processing computed tomographic data, DICOM files and generating Shepp-Logan phantoms.

SSH clone
```bash
git clone git@github.com:kulvait/KCT_denpy.git
```

## Instalation

To install the package, execute the following command

```bash
pip install git+https://github.com/kulvait/KCT_denpy.git
```

For an editable local install from the git directory, use the following command

```bash
git clone https://github.com/kulvait/KCT_denpy.git
cd KCT_denpy
pip install --user --upgrade -e .
```


## Upgrading the Package
To update the package, use

```bash
pip install --upgrade git+https://github.com/kulvait/KCT_denpy.git
```

For a local upgrade from the git directory:

```bash
pip install --user --upgrade .
```

## Modules Overview

### DEN
The denpy.DEN module provides tools to manipulate DEN files in Python. The DEN format includes a short header specifying the dimensionality of the data, making it particularly convenient for computed tomography applications.

**Typical import of the module**
```python
from denpy import DEN
```

**Getting information about the DEN file**
- `info = DEN.readHeader(file_path)`: Reads the header of a DEN file and returns a dictionary with the header information.

**Important fields of the info dictionary**

- `info["dimspec"]`: Dimensions of the data (e.g., (512, 256, 128) for a 3D dataset, which corresponds to 128 frames of xdim x ydim = 512 x 256). 
- `info["dimcount"]`: Number of dimensions in the data `info["dimcount"] = len(info["dimspec"])`.
- `info["shape"]`: Shape of the array that can be read with `DEN.getNumpyArray`, which is flipped compared to `info["dimspec"]` because of the way NumPy indexes arrays.

**Reading and writing entire files**

- `x = DEN.getNumpyArray(file_path)`: Reads the entire DEN file and returns it as a NumPy array.
- `DEN.writeEmptyDEN(file_path, [dimx, dimy, dimz], force=True)`: Writes an empty DEN file with the specified dimensions.
- `DEN.storeNdarrayAsDEN(file_path, x, force=True)`: Writes a NumPy array to a DEN file.

Note that flipped order of the Numpy array dimensions allows user to use f = x[0] to read the first frame in the file, keep in mind that first index to the Numpy array is the slowest changing index, so the last index in `info["dimspec"]`.

**Reading individual frames**

- `DEN.getFrame(file_path, index)`: Reads a single frame (slice) from a DEN file and returns it as a NumPy array.
- `DEN.writeFrame(file_path, index, frame, force=True)`: Writes a single frame (slice) to a DEN file.

### ZAR, ZARZS, ZARPYRAMID
Zarr helpers: `ZAR.get_compressor` returns zarr v2/v3 codecs by name (Blosc with bit shuffle, zstd,
lz4, gzip and the image codecs htj2k, jpeg2k, jpegxl, avif, ...); `ZARZS.ZipStore` is a drop-in for
`zarr.storage.ZipStore` that threads write and read in parallel (as fast as a directory store), with
deletes and safe appends (`ZARZS.recover`); `ZARPYRAMID` writes multiscale pyramids with OME-NGFF
metadata.

### TIFF
Stacks of TIFF files to zarr or DEN, described by a JSON specification: tiff stacks (which files,
frame shape, dtype, known missing frames), virtual stacks computed with lambdas or reductions
(e.g. a dark median and dark corrected projections) and the arrays to write. Everything is validated
before writing; frames are read, computed, compressed and written in parallel. The reference of the
format is [docs/TIFFSTACK_SPEC.md](docs/TIFFSTACK_SPEC.md).

```bash
pip install "denpy[zarr]"                     # zarr, tifffile, imagecodecs
convertTiffStack scan.json --dry-run          # validate and print what would be written
convertTiffStack scan.json -j 32 --report report.json
convertTiffStack --files img_*.tif --zarr img.zar --compression htj2k --clevel 0
```

```python
from denpy import TIFF
report = TIFF.execute(TIFF.validate("scan.json"))
```

### PETRA (legacy)
Parsing of PETRA III scan HDF5 files (`scanDataset`, `getExperimentInfo`, ...). Kept for the existing
callers; the object oriented implementation is `lautools.PETRA` (`from lautools import PETRA`, same
module functions), described in the lautools API manual.

### COR
The denpy.COR module is used to detect the center of rotation in tomographic data from synchrotron sources.

### DICOM
The denpy.DICOM module helps in reading DICOM files and is mainly intended for parsing series of data from CT perfusion acquisitions.

### PHANTOM
The denpy.PHANTOM module generates various forms of the Shepp-Logan phantom. It aims to provide a well-defined version of this phantom in both 2D and 3D, addressing inconsistencies in previous definitions.

### UTILS
The denpy.UTILS module provides utilities for parsing simple parameter files.

## Example usage

The following code demonstrates the basic usage of the `denpy` package for reading a stack of 3D projection data and converting them to TIFF images.

```python
from denpy import DEN # Import the DEN module
from PIL import Image
import numpy as np
import os

def write_slice(img, file_name, force=False):
    if os.path.exists(file_name) and not force:
        raise IOError(f"File {file_name} already exists!")
    im = Image.fromarray(img, mode='F')  # float32
    im.save(file_name, "TIFF")

def convert_den_to_tiff(input_den, output_dir, force=False, suffix=""):
    header = DEN.readHeader(input_den)
    if header["dimcount"] < 2:
        print(f"File {input_den} is {header['dimcount']}-dimensional, too few dimensions to convert to TIFF.")
        return

    xdim, ydim = header["dimspec"][:2]
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    if header["dimcount"] == 2:
        write_slice(DEN.getFrame(input_den, 0), f"{output_dir}/{suffix}.tif", force)
    else:
        frame_count = np.prod(header["dimspec"][2:])
        significant_digits = int(np.log10(frame_count)) + 1
        format_string = f"{{:0{significant_digits}d}}"

        for index in np.ndindex(header["dimspec"][2:]):
            index_str = "_".join(format_string.format(e) for e in index)
            write_slice(DEN.getFrame(input_den, index), f"{output_dir}/{index_str}{suffix}.tif", force)

# Example usage
input_den = 'path/to/input.den'
output_dir = 'path/to/output/dir'
convert_den_to_tiff(input_den, output_dir, force=True, suffix='_projection')
```

### Explanation
1. `write_slice` function: Saves an image slice as a TIFF file, optionally overwriting existing files.
2. `convert_den_to_tiff` function:
 - Reads the header of the DEN file to get dimensional information.
 - Checks if the data has at least two dimensions.
 - Creates the output directory if it doesn't exist.
 - Writes each frame (slice) to a TIFF file, formatting the file names according to the slice index.

### To use the example code:
Replace 'path/to/input.den' and 'path/to/output/dir' with the actual paths to your DEN file and desired output directory. The `force=True` argument allows overwriting existing files, and suffix='__projection' adds a suffix to the output file names. Adjust these as needed for your use case.


## Example scripts

In the `scripts` directory, you can find example scripts and useful tools for working with the denpy package. These files are not part of the main package but can be helpful references.

## Documentation

The manual [manual/denpy_manual.pdf](manual/denpy_manual.pdf) describes DEN files, the zarr helpers
and the TIFF stack conversion. Build it with `tectonic manual/denpy_manual.tex`.

## Tests

```bash
python -m pytest denpy/
```

## Licensing

Unless otherwise specified in the source files, this project is licensed under the GNU General Public License v3.0.

Copyright (C) 2018-2026 Vojtěch Kulvait

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, version 3 of the License.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.

# Donations

If you find this software useful, you can support its development through a donation.

[![Thank you](https://img.shields.io/badge/donate-$15-blue.svg)](https://kulvait.github.io/donate/?amount=15&currency=USD)
