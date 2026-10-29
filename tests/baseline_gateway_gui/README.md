# gateway_gui refactor verification baseline

Snapshot of `pyMT.GatewayGUI.gateway_main` behaviour, captured on the
`cleanup/gateway_gui` branch **before** refactoring.

`gateway_main` is a thin PyQt launcher (a dataset table, save / load of a
`.pymt` project file, and three buttons that spawn the other GUIs). The
spawning itself is covered by the individual GUI baselines and the manual
smoke-test checklist (`ModelWindow` segfaults head-less on this machine).
This harness freezes the module's *pure* logic, which runs fine under
`QT_QPA_PLATFORM=offscreen`. `tests/gateway_gui_baseline.py` builds a
real `NewProject` with the open-project dialog stubbed and freezes:

* construction state — table row / column count, header labels,
  `empty_table`, `get_locs_from()`, the `file_to_index` / `index_to_file`
  maps,
* `sort_files` extension classification for real ModEM files plus a
  spread of extension-only and unknown-extension inputs,
* the dataset-table round trip — `init_dataset_table` then
  `dataset_dict()` / `active_datasets()` (per selected row) /
  `dataset_index()` / `empty_table`,
* the text `write_project` emits for a known two-dataset table, and the
  `FileInputParser.read_pystart` parse of that same file.

## Notes captured by the baseline (pre-existing, not changed here)

* `sort_files` classifies `.dat` data-vs-response by looking for
  `"response"` in the first line; pyMT-written ModEM response files do
  not contain that word, so `swz_test_modem_resp.dat` is filed as
  `data`.
* `write_project` writes `resolution <path>` lines, but
  `FileInputParser.read_pystart` does not list `resolution` in its
  accepted keywords, so it prints `Unrecognized keyword: resolution` and
  drops the value on load.

## Files

| File | What it is |
|---|---|
| `gateway_gui_baseline.pkl` | pickled `{label: value}` snapshot (23 entries) |
| `gateway_gui_baseline_summary.txt` | human-readable shapes / crc32s |
| `scratch/` | the written `.pymt` project file (git-ignored) |

## Regenerate / compare

```
QT_QPA_PLATFORM=offscreen python tests/gateway_gui_baseline.py save
QT_QPA_PLATFORM=offscreen python tests/gateway_gui_baseline.py compare
```
