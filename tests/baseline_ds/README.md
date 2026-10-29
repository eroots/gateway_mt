# data_structures refactor verification baseline

Snapshot of `pyMT.data_structures` public behaviour, captured on the
`cleanup/data_structures` branch **before** refactoring, for before/after
comparison.

## Files

| File | What it is |
|---|---|
| `ds_baseline.pkl` | pickled `{label: value}` snapshot, ~1260 entries |
| `ds_baseline_summary.txt` | human-readable shapes / sample values |
| `roundtrip_*.dat` | scratch files from the write -> re-read round-trip checks |

## Regenerate / compare

From the repo root, with the `pymt` conda env:

```
python tests/data_structures_baseline.py save      # overwrite the snapshot
python tests/data_structures_baseline.py compare   # re-run and diff
```

`compare` exits non-zero and lists every entry that changed shape, value
(outside `numpy.allclose` tolerance), or exception status.

## Inputs

* `test_data/claude_test_data/edis/allall.lst` (85 EDIs) -> `RawData`
* `test_data/claude_test_data/swz_test_modem.dat` -> `Data` / `Dataset` (ModEM)
* `test_data/claude_test_data/swz_test_gofem.dat.gdat` -> `Data` (GoFEM)
* `test_data/testfile_WS.data` -> `Data` (WSINV3DMT)
* `test_data/testfile_MARE2DEM.emdata` -> `Data` (MARE2DEM)
* `test_data/testfile_ModEM.model` -> `Model`

## Coverage

* `RawData`: read, `NS`, `site_names`, `locations`, `origin`, `master`/
  `narrow_period_list`, `get_locs` at several azimuths, `average_rho`,
  per-site data/error/errmap dicts, phase-tensor parameters, `get_data`.
* `Data` (4 formats): `NP`/`NS`/`NR`, `periods`, `components`, `inv_type`,
  `used_components`, `has_flagged_data`, `get_locs`, `average_rho`, per-site
  dicts; on copies: `rotate_sites`, `apply_error_floor`, `calculate_error_floor`,
  `detect_outliers`, `Site.get_data`.
* `Model`: `nx/ny/nz`, `dx/dy/dz`, `xCS/yCS/zCS`, `center`, `vals` stats,
  `cell_centers`, `is_half_space`, `conductance`, `depth_weighted_rho`,
  `elevation`; on copies: `dz_insert`, `dx_delete`, `to_local`.
* `Dataset`: `azimuth`, `site_names`, `calculate_RMS` (against a
  deterministically perturbed pseudo-response).
* `PhaseTensor.__sub__` for all six `difference_definition` values.
* write -> re-read round-trip for the ModEM `Data`.

## Known / expected

* `mare/average_rho` records a `KeyError` -- the MARE2DEM sample data has no
  `ZXYR` component, so `average_rho` (which needs the full impedance tensor)
  legitimately fails. `compare` matches EXCEPTION == EXCEPTION.
* No model response file for `swz_test_modem.dat` was available, so
  `calculate_RMS` is baselined against a perturbed copy of the data.
