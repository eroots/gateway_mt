# data_plot refactor verification baseline

Snapshot of `gateway_mt.DataGUI.data_plot` behaviour, captured on the
`cleanup/data_plot` branch **before** refactoring.

`data_plot` is a PyQt GUI, so `tests/data_plot_baseline.py` runs it
head-less (`QT_QPA_PLATFORM=offscreen`): it builds a real `DataMain`
from `tests/data/synthetic/synthetic.pymt` (the synthetic list / data /
response / model), drives a fixed script of ~30 menu / toggle /
sort / flag / coordinate / write handlers with every modal dialog
monkey-patched to a canned answer, and freezes:

* `DataPlotManager` config (`components`, `scale`, `errors`, `toggles`,
  `wrap_phase`, `outlier_thresh`, tiling, per-axis line xy-data and
  y-limits) before and after the script,
* `MapView` config (all `*_cax`, colour map, `coordinate_system`,
  marker / PT / induction params, `spatial_units`, `interpolant`) and the
  three `site_locations` arrays (rows sorted — `generic` order is a set
  difference),
* `Dataset` state — NP/NS/NR, `inv_type`, azimuth, component list, a few
  sample sites' periods / data / `used_error`, station RMS,
* the text of every file the GUI writes (list + ModEM + ModEM-removed +
  WSINV3DMT + GoFEM); files over 100 kB are stored as
  `crc32 + line count + head/tail` rather than in full,
* one entry per driven step recording whether it raised.

## Files

| File | What it is |
|---|---|
| `data_plot_baseline.pkl` | pickled `{label: value}` snapshot (385 entries) |
| `data_plot_baseline_summary.txt` | human-readable shapes / crc32s |
| `scratch/` | files the write handlers produced (git-ignored) |

## Regenerate / compare

From the repo root, with the `pymt` conda env:

```
QT_QPA_PLATFORM=offscreen python tests/data_plot_baseline.py save
QT_QPA_PLATFORM=offscreen python tests/data_plot_baseline.py compare
```

`compare` exits non-zero and lists every entry whose value moved outside
tolerance. A full run takes about a minute (loads the raw EDIs and redraws the
map several times).

## History

The original snapshot (taken on `cleanup/data_plot` before the gplot
fixes were merged) recorded `step/map_planview` as an `AttributeError`
from gplot's misspelled `use_colorbar`. With gplot fixed, the current
snapshot has no exception entries.

## Not a substitute for a manual pass

This harness proves the handlers still *run* and the written files stay
byte-identical; it cannot see the rendered figures or true click / drag
interaction. Pair it with the manual smoke-test checklist in the module
changelog.
