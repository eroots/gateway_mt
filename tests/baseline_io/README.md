# IO refactor verification baseline

Snapshot of `pyMT.IO` reader/writer behaviour, captured on the `cleanup/IO`
branch **before** refactoring.

## Files

| File | What it is |
|---|---|
| `io_baseline.pkl` | pickled `{label: value}` snapshot (~505 entries) |
| `io_baseline_summary.txt` | human-readable shapes / crc32 of written text |
| `scratch/` | files written by the round-trip / VTK / CSV checks (git-ignored) |

## Regenerate / compare

From the repo root, with the `pymt` conda env:

```
python tests/io_baseline.py save      # overwrite the snapshot
python tests/io_baseline.py compare   # re-run and diff
```

`compare` exits non-zero and lists every entry whose parsed value or
written text changed.

## Coverage

* `read_data` for **modem / gofem / wsinv3dmt / mare2dem** -- `other_info`
  plus per-site `data` / `errors` / `periods` / `locations` for the sample
  sites.
* `read_raw_data` from the 30 synthetic EDIs.
* `read_model` (`synthetic.model`), `read_data_header` (modem/ws/mare),
  `read_sites`, `read_freqset`, `get_components`.
* `write_data` round-trips (write -> re-read -> compare parsed data, plus a
  normalised copy of the written text) for all four data formats.
* `write_model` round-trip; `write_covariance` -> `read_covariance` round-trip.
* `model_to_vtk`, `sites_to_vtk`, `write_phase_tensors`,
  `write_induction_arrows`, `write_locations` (csv), `write_edi`,
  `write_response`, `write_list` -- written text frozen (dates / scratch
  paths blanked).

Written text is compared verbatim after blanking timestamps and the
scratch directory path.

## Expected diffs after the bug-fix pass

`compare` will report exactly one intended change:

* `roundtrip/modem/written_text` line 1: `UTM Zone: {}` -> `UTM Zone: Undefined`
  -- the ModEM writer's header `.format()` was binding to the wrong
  concatenated string, so files were written with a literal `{}`. Now the
  actual `data.UTM_zone` is substituted.

Everything else stays byte-identical.

## Not covered (no fixtures)

* `read_covariance` exception-block path (needs a ModEM covariance file with
  air/ocean masks) -- also where the `model.nx` NameError lives.
* EM3DANI / MT3DANI models, `.adat` data, Occam 2D.
* `write_locations` shapefile / kml branches.
