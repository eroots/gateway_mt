# utils refactor verification baseline

Snapshot of `pyMT.utils` public behaviour, captured on the `cleanup/utils`
branch **before** any refactoring, for before/after comparison.

## Files

| File | What it is |
|---|---|
| `utils_baseline.pkl` | pickled `{label: value}` mapping, 505 entries — the snapshot |
| `utils_baseline_summary.txt` | human-readable shapes / sample values for every entry |
| `geotools_filter.png` | `geotools_filter` output on a synthetic spiky curve |
| `compute_bost1D.png` | `compute_bost1D` depth curve for one EDI site |

## Regenerate / compare

From the repo root, with the `pymt` conda env:

```
python tests/utils_baseline.py save      # overwrite the snapshot
python tests/utils_baseline.py compare   # re-run and diff against the snapshot
```

`compare` exits non-zero and lists every entry that changed shape, value
(outside `numpy.allclose` tolerance), or exception status.

## Inputs

* `test_data/claude_test_data/edis/allall.lst` (85 EDI files) -> `RawData`
* `test_data/claude_test_data/swz_test_modem.dat` -> ModEM `Data`
* `test_data/claude_test_data/swz_test_gofem.dat.gdat` -> GoFEM `Data`
* `test_data/testfile_ModEM.model` -> `Model` (for `regrid_model`)
* frozen synthetic arrays defined in `_frozen_synthetic_inputs()`

## Coverage notes

* Every public function in `utils.py` is exercised.
* The snapshot was captured **before** the refactor, so it still contains
  entries for the functions that were then removed as dead code
  (`rms`, `strike_direction`, `enforce_output`, `row_or_col`,
  `remove_bad_data`, `normalize_resolution`).  `compare` reports those as
  "Missing now (were in baseline)" — that is expected.
* `dms_to_dd` is the one intentional value change: it is now an alias of
  `dms2dd`, so `compare` reports a `max abs diff` for it. Also expected.
* `calculate_misfit` is fed a real data site vs. a deterministically perturbed
  copy of itself as the "response".
* `crosscheck/*` runs real EDI lat/longs through `project` and `geo2utm` and
  confirms they agree with each other and with `RawData`'s own projection.
