# mesh_designer refactor verification baseline

Snapshot of `pyMT.ModelGUI.mesh_designer` behaviour, captured on the
`cleanup/mesh_designer` branch **before** refactoring.

Unlike `model_viewer`, `mesh_designer` is a pure-matplotlib PyQt GUI and
runs head-less (`QT_QPA_PLATFORM=offscreen` + Agg).
`tests/mesh_designer_baseline.py` builds a real `model_viewer_2d` from
`synthetic.model` + `synthetic_modem.dat`, drives a fixed script of
handlers with the modal dialogs monkey-patched, and freezes:

* the plan-view `pcolormesh` array / clim and the overlaid line data,
  axis limits / aspect, the NX/NY/NZ labels,
* the model cell-edge arrays (`dx` / `dy` / `dz`), `(nx, ny, nz)`, and
  `vals` statistics, before and after the script,
* the text of the written ModEM model + covariance,
* one entry per driven step recording whether it raised.

Driven steps: z-slice change, colour map, mesh-line colour, background
resistivity, Gaussian smooth, add / remove padding, insert Y / X mesh
lines by synthetic click, Ctrl-click delete, min/max depth edits, lock
aspect ratio, write ModEM.

## Not driven

`regenerate_mesh` (`Model.generate_mesh`) hung on the original 134-site
fixture, and `generate_depths_decades`
(`Model.generate_zmesh`) needs a populated decade list. Both are heavy
`data_structures` operations, out of scope for the mesh_designer
docstring pass.

## Files

| File | What it is |
|---|---|
| `mesh_designer_baseline.pkl` | pickled `{label: value}` snapshot (46 entries) |
| `mesh_designer_baseline_summary.txt` | human-readable shapes / crc32s |
| `scratch/` | the written model / covariance (git-ignored) |

## Regenerate / compare

```
QT_QPA_PLATFORM=offscreen python tests/mesh_designer_baseline.py save
QT_QPA_PLATFORM=offscreen python tests/mesh_designer_baseline.py compare
```
