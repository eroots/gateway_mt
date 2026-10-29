# pyMT/scripts

Command-line tools installed with `pyMT` (via `setup.py`). Nothing in
`pyMT` imports from this folder; only `j2ws3d.py` imports *out* of it
(to `pyMT.ModelGUI.mesh_designer`).

Personal analysis and figure scripts that used to live here are now in
[`examples/`](../../examples/) at the top of the repository.

## Shipped / CLI tools

Installed by `setup.py` (`scripts=[...]` and/or `console_scripts`). These
take command-line arguments or a small parameter block and are meant to
be run as-is.

| Script | Entry | Purpose |
|---|---|---|
| `f2l.py` | `console_scripts: f2l` | Build a `.lst` site-list file from a folder of EDI / `.dat` files. |
| `to_vtk.py` | `console_scripts: to_vtk` + `scripts` | Write a `Data` / `Model` (+ locations) to VTK for ParaView. |
| `plot_edi.py` | `console_scripts: plot_edi` | Plot apparent resistivity / phase / tipper curves from EDI files. |
| `ws2modem.py` | `scripts` | Convert a WSINV3DMT data file to ModEM format. |
| `plot_occam.py` | `scripts` | Plot Occam-2D data and response. |
| `autogen_report.py` | `scripts` | Scan the working directory for model files and build a PDF misfit report. |
| `j2ws3d.py` | `scripts` | Convert j-format `.dat` data to WSINV3DMT, and (with `-i`) launch the mesh designer. |
| `plot_dataset.py` | `scripts` | Plot a whole dataset (data / response / raw) to a multi-page PDF via `gplot.DataPlotManager`. Paths are edited in the file, not passed on the command line. |
| `plot_transect.py` | `scripts` | Small transect-location plot from a shapefile. Paths hard-coded in the file. |
