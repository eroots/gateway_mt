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

All from `tests/data/synthetic/`:

* `edis/synthetic.lst` (30 EDIs) -> `RawData`
* `synthetic_modem.dat` -> `Data` / `Dataset` (ModEM)
* `synthetic_gofem.gdat` -> `Data` (GoFEM)
* `synthetic_ws.data` -> `Data` (WSINV3DMT)
* `synthetic_mare2dem.emdata` -> `Data` (MARE2DEM)
* `synthetic.model` -> `Model`

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

* `mare/load` records a `ValueError`: pyMT's own MARE2DEM writer puts
  single spaces between the UTM zone and origin values in the header,
  but the reader splits that line on two or more spaces, so pyMT cannot
  re-read a MARE2DEM file it wrote with a real origin (known bug, not
  fixed here). `compare` matches EXCEPTION == EXCEPTION.
* `calculate_RMS` is baselined against a deterministically perturbed copy
  of the data.
