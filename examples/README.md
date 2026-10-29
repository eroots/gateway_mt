# examples

Scripts that use `pyMT` for specific analyses and figures. They are not
installed with the package -- copy one, edit the paths and options at the
top, and run it. Most were written for a particular survey, so expect to
adapt them. See each file's module docstring for details.

The installed command-line tools are in [`pyMT/scripts/`](../pyMT/scripts/).

## Analysis and figure scripts

"Edit the paths at the top and run" scripts -- a research record, not
general tools. Most have absolute paths baked in and top-level code with
no `main()`. Several need packages **not** in `pyMT`'s
`install_requires`.

| Script | Extra deps | What it does |
|---|---|---|
| `1D_response.py` | -- | 1-D layered-earth MT forward response (functions only). |
| `apply_distortion.py` | -- | Apply a frequency-dependent distortion tensor to a data set. |
| `apply_distortion_to_EDI.py` | `mtpy` | Same, writing distorted EDIs. |
| `archies_law.py` | -- | Archie's-law / Hashin-Shtrikman mixing curves. |
| `average_1D_model.py` | -- | Average many 1-D models and compare their responses. |
| `bath2model.py` | `rasterio`, `cmocean` | Add bathymetry / topography to a model from an SRTM raster. |
| `calc_aniso_rho.py` | -- | Estimate anisotropic resistivities from an isotropic slice (Eisel & Haak 1999). |
| `calculate_distortion_matrix.py` | -- | Recover the per-site distortion tensor between two data sets. |
| `castro_ellipse_test.py` | -- | Scratch sketch of a rotated ellipse with axes / arrows (does not use `pyMT`). |
| `conductance_to_geotiff.py` | `rasterio` | Depth-integrated conductance of a model to GeoTIFF. |
| `copy_data.py` | -- | Add raw-data tipper to a GoFEM data file, taking tipper errors from two other data files. |
| `edi_fix_rotation.py` | `mtpy` | Interpolate / de-rotate a folder of EDIs. |
| `edi_spec_to_tf.py` | `mtpy` (old) | Re-write spectra EDIs as transfer-function EDIs. |
| `generate_colourbar.py` | -- | Dump every `e_colours` colour map to a PNG. |
| `gofem_min_dists.py` | -- | List receiver pairs closer than a threshold. |
| `image_to_vtk.py` | `PIL` | Georeferenced raster (`.tif` + world file) to a VTK surface. |
| `interp_plan_rms.py` | `naturalneighbor` | Interpolated plan-view RMS map / scatter. |
| `interp_rhopha.py` | `naturalneighbor` | Interpolated rho / phase plan sections and differences. |
| `lst2kml.py` | `simplekml`, `pyshp` | Site list to KML / shapefile. |
| `model_to_geotiff.py` | `rasterio` | Model depth slices / conductance to GeoTIFF. |
| `np_to_geotiff.py` | `rasterio` | Grid lon / lat / z points to a GeoTIFF (has a real `argparse` CLI). |
| `paint_model.py` | -- | Hand-edit resistivity regions of a model. |
| `paint_anisotropic_model.py` | -- | Same for an anisotropic (EM3DANI) model. |
| `phase_tensor_map.py` | -- | Plan-view phase-tensor ellipse map (via `gplot.MapView`). |
| `plan_pseudosection.py` | `naturalneighbor` | Interpolated plan pseudosection. |
| `plot_bathymetry.py` | `basemap` | Download + plot NOAA bathymetry for a lat/long box. |
| `plot_dataset_rhophatip.py` | -- | Rho / phase / tipper dataset plots (self-contained, does not use `gplot`). |
| `plot_induction_map.py` | -- | Plan-view induction-arrow map (via `gplot.MapView`). |
| `plot_misfit_curves.py` | -- | Per-component misfit-curve PDFs (function `generate_misfit_plots`). |
| `plot_slice_example.py` | -- | Worked example of interpolating an arbitrary slice through a model. |
| `pseudosection_gplot.py` | -- | Batch pseudosections via `gplot`. |
| `pt2csv.py` | `pandas` | Phase-tensor parameters for every site/period to CSV. |
| `pt_pseudosection.py` | -- | Phase-tensor ellipse pseudosection along a profile. |
| `read_csr.py` | -- | Read a MARE2DEM CSR penalty-matrix file (does not use `pyMT`). |
| `read_shape_plot.py` | `cartopy`, `pyshp` | Overlay a shapefile on station locations. |
| `rewrite_edis.py` | -- | Re-write EDI files from pyMT `Site` objects. |
| `rms_histogram.py` | -- | RMS-misfit histograms / difference histograms. |
| `rms_period_barchart.py` | -- | Per-period RMS-misfit curves comparing several responses. |
| `scatter_bostick.py` | -- | 3-D scatter of Niblett-Bostick depth vs tipper. |
| `shapefile_points_to_vtk.py` | `geopandas`, `pyvista` | Shapefile geometry to VTK. |
| `test_bosts.py` | -- | Development scratch: a local re-implementation of `compute_bost1D` with test plots. |
| `to_vtk_point-cloud.py` | -- | Like `to_vtk.py` but writes a point cloud. |
| `write_induction_arrows.py` | -- | Write an induction-arrow table (template: edit three variables). |
| `write_phase_tensors.py` | -- | Write a phase-tensor ellipse table. |
| `ztem_to_edi.py` | -- | Convert ZTEM Geosoft `.grd` grids to EDI files. |

## Overlapping / superseded

Groups of scripts that do nearly the same thing, or whose job the GUIs
now cover. Kept for reference; each carries a header note pointing at
its replacement.

* **Crooked-transect model slices**: `plot_model_slice.py`,
  `integrated_crooked_slice.py`, `crooked_model_slice_commented.py` --
  three overlapping versions; `model_viewer`'s transect tab now does
  this interactively.
* **Pseudosections**: `plan_pseudosection.py`, `pseudosection_gplot.py`,
  `pt_pseudosection.py`, `interp_rhopha.py`, `interp_plan_rms.py` --
  overlap with `gplot.MapView.plan_pseudosection` and `data_plot`'s map
  tab.
* `pyMT/scripts/plot_dataset.py` vs `plot_dataset_rhophatip.py` -- two
  dataset plotters.
* `conductance_to_geotiff.py` vs `model_to_geotiff.py` -- both note the
  calculation is now a `Model` method.
* `pyMT/scripts/to_vtk.py` vs `to_vtk_point-cloud.py`.

## Candidates to promote

If these are turned into real CLIs or folded into `data_plot` later,
they are the natural starting points (already parameterised, minimal
personal state): `plot_misfit_curves.py`, `write_induction_arrows.py`,
`write_phase_tensors.py`, `ztem_to_edi.py`, `plot_dataset_rhophatip.py`.
