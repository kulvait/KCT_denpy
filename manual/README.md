# denpy documentation

- `denpy_manual.pdf` – DEN files, zarr helpers (`ZAR`, the parallel zip store `ZARZS`,
  `ZARPYRAMID`), the TIFF stack conversion (`TIFF`, `convertTiffStack`) and the legacy `PETRA`.
- `../docs/TIFFSTACK_SPEC.md` – the reference of the TIFF stack specification.

Build the PDF with [tectonic](https://tectonic-typesetting.github.io/):

```bash
tectonic denpy_manual.tex
```

Scan metadata and the conversion of PETRA III scans are in lautools (`manual/lautools_api_manual.pdf`
of lautools).
