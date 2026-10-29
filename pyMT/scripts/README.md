# pyMT/scripts

Command-line tools installed with `pyMT` (via `setup.py`). Nothing in
`pyMT` imports from this folder.

Personal analysis and figure scripts that used to live here are now in
[`examples/`](../../examples/) at the top of the repository.

## Shipped / CLI tools

Installed by `setup.py` (`scripts=[...]` and/or `console_scripts`). These
take command-line arguments or a small parameter block and are meant to
be run as-is.

| Script | Entry | Purpose |
|---|---|---|
| `f2l.py` | `console_scripts: f2l` | Write a site list (`all_sites.lst`, or `<FILTER>_sites.lst`) of the EDI / `.dat` files in a folder. Run `f2l -h` for options. |
| `to_vtk.py` | `console_scripts: to_vtk` + `scripts` | Write a model and / or site locations to VTK for ParaView (ModEM / WSINV3DMT-style inputs). Interactive with no arguments, or `to_vtk -h` for the command-line options. |
| `plot_edi.py` | `console_scripts: plot_edi` | Plot apparent resistivity / phase / tipper curves from EDI files: shown on screen, or saved as PNG / PDF (`plot_edi all`, `--save`). See `plot_edi -h`. |
| `autogen_report.py` | `scripts` | Misfit report for an inversion (any data format pyMT reads): `summary.txt` with RMS by component and worst sites, plus `misfit_curves.pdf` (rho / phase or impedance, with tipper; four sites per page). See `autogen_report -h`. |
| `plot_dataset.py` | `scripts` | Plot a whole dataset (data / response / raw) to a multi-page PDF via `gplot.DataPlotManager`. Paths are edited in the file, not passed on the command line. |
| `plot_transect.py` | `scripts` | Small transect-location plot from a shapefile. Paths hard-coded in the file. |
