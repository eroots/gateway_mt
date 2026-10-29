# Tests

Run the suite from the repository root:

```
pip install -e .[test]
pytest
```

The GUI tests use Qt's offscreen platform, so no display is needed. A full
run takes about half a minute.

| File | What it checks |
|---|---|
| `test_physics.py` | Apparent resistivity, phase and phase tensors against known results (a uniform half-space, an independent 1-D recursion); rotation and its invariants; coordinate conversions against pyproj. |
| `test_io.py` | Reading EDIs and site lists; write / read round trips for each data format and for models; the same data in different formats agrees; the response file holds the analytic responses it was built from. |
| `test_data_structures.py` | `Dataset` / `Data` / `Model` behaviour on the synthetic survey: RMS ~1, flagging of periods outside a station's band, resampling raw data, rotation, adding / removing sites, mesh edits. |
| `test_plotting.py` | Sounding-curve plots show the computed values; the map view draws one arrow / ellipse per station; pseudosections. |
| `test_commands.py` | `f2l`, `plot_edi`, `plot_dataset`, `autogen_report` and `to_vtk` run on the synthetic data. |
| `test_gui.py` | The gateway window, Data Plot (incl. writing data), its map view and Mesh Designer, offscreen; Model Viewer's model-to-VTK conversion. |
| `test_synthetic_data.py` | Regenerating the synthetic dataset reproduces the committed files -- an end-to-end check of the readers and writers. |

`data/synthetic/` holds the test data and `data/make_synthetic.py` the
script that generates it (see the README there).

## Known bugs

Tests marked `xfail` document known bugs: they are expected to fail until
the bug is fixed, and the reason says what is wrong. Once a fix lands the
test passes, pytest reports it as an unexpected pass (`XPASS(strict)`),
and the marker should be removed.
