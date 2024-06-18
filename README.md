# denpy

denpy is a Python package designed for manipulating the less well-known DEN (Dennerlein) raw data format, commonly used in computed tomography. This package also includes modules for working with other data formats and processing tasks related to tomography.

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

## Licensing

Unless otherwise specified in the source files, this project is licensed under the GNU General Public License v3.0.

Copyright (C) 2018-2024 Vojtěch Kulvait

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
