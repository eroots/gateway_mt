# model_viewer refactor verification baseline

Snapshot of `pyMT.ModelGUI.model_viewer` behaviour, captured on the
`cleanup/model_gui` branch **before** refactoring.

## Why only `model_to_rectgrid`

`ModelWindow` is a PyVista (`pyvistaqt.QtInteractor`) GUI. Constructing
it head-less on this machine segfaults inside VTK ("failed to get valid
pixel format" -- no usable OpenGL 3.2 context), so an automated GUI
harness isn't possible here. Per CLAUDE.md, GUIs are verified with a
**manual smoke-test checklist** instead (see the module changelog).

What *is* automatable is the module's one reusable, render-free entry
point: `model_to_rectgrid(model, resolution=None, rho_axis='rho_x')`.
This harness freezes its output for `tests/data/synthetic/synthetic.model`:

* the `pyvista.RectilinearGrid` `dimensions`, `n_cells`, `bounds`,
* every `cell_data` array (`Resistivity`, `Resolution`, and
  `RawResolution` when a resolution volume is passed),

with and without a synthetic resolution volume (the model's own cell
sizes as a positive scalar field, to exercise the normalisation branch).

`rho_y` / `rho_z` and the anisotropic-ratio `rho_axis` options aren't
covered -- the sample model is isotropic and there is no anisotropic
fixture. (`model_to_rectgrid` raises on `'rho_x/rho_y'` for an isotropic
model; the GUI never reaches that path because `change_rho_axis` is only
connected for anisotropic models.)

## Files

| File | What it is |
|---|---|
| `model_viewer_baseline.pkl` | pickled `{label: value}` snapshot (16 entries) |
| `model_viewer_baseline_summary.txt` | human-readable shapes / sample values |

## Regenerate / compare

```
python tests/model_viewer_baseline.py save
python tests/model_viewer_baseline.py compare
```

`compare` exits non-zero and lists any entry that moved outside
`numpy.allclose` tolerance.
