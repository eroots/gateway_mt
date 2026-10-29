# gplot refactor verification baseline

Snapshot of `pyMT.gplot` rendering behaviour, captured on the
`cleanup/gplot` branch **before** refactoring.

`gplot` returns matplotlib artists, not numbers, so the harness drives
every public `DataPlotManager` / `MapView` method against the sample data
and freezes the *artist state* each one produces:

* line `xy`-data and colour,
* collection offsets / colour arrays / clims / quiver `UVC` (all row-sorted
  so set-ordering in site lists can't cause spurious diffs),
* patch vertices and facecolours,
* image arrays / clims / extents,
* axis title, x/y labels, limits, scale, aspect,
* `format_coord` readout strings at a fixed cursor point,
* the pure-numeric helpers (`pt_fill_limits`, `generate_*`, `get_label`,
  `interpolate`, `get_tensor_params`, `gen_tiling`, `units`,
  `get_designated_colour`, …).

## Files

| File | What it is |
|---|---|
| `gplot_baseline.pkl` | pickled `{label: value}` snapshot (~8500 entries) |
| `gplot_baseline_summary.txt` | human-readable shapes / sample values |

## Regenerate / compare

From the repo root, with the `pymt` conda env:

```
python tests/gplot_baseline.py save      # overwrite the snapshot
python tests/gplot_baseline.py compare   # re-run and diff
```

`compare` exits non-zero and lists every entry whose value moved outside
tolerance (`1e-9` for pure arithmetic, `1e-6` for anything routed through
pyproj / interpolation / geometry).

## Inputs

* `test_data/claude_test_data/edis/allall.lst` (+ the SWZ EDIs) -- raw data
* `test_data/claude_test_data/swz_test_modem.dat` -- inverted data
* `test_data/claude_test_data/swz_test_modem_resp.dat` -- synthetic response
  (data + 10 % noise)
* `test_data/testfile_ModEM.model` -- model for the plan / slice plots
  (does not spatially overlap the SWZ survey; used only for reproducible
  coverage of the model-plotting code paths)

`interpolant` is pinned to `linear` and `numpy`'s RNG is seeded so the run
is deterministic; colourbars are left off (they mutate figure-level state
that isn't worth snapshotting).

## Expected diff after the pass-5 bug-fix commit

`compare` reports exactly one intended change:

| Entry | Baseline | After pass 5 |
|---|---|---|
| `dpm/redraw_single/no_1d_key` | `EXCEPTION: KeyError: '1d'` | runs cleanly | `redraw_single_axis` no longer indexes a missing `self.sites['1d']` key. |

Everything else stays byte-identical (the other pass-5 fixes touch paths
that are either unreachable from this harness or produce the same
output).

`plot_plan_view` originally read `self.use_colorbar` (missing the `u`);
the harness sets both spellings to `False` so the model plots run either
way, which is why fixing the typo does not move the snapshot.
