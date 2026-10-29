"""Regenerating the synthetic dataset reproduces the committed files exactly.

make_synthetic.py reads its EDIs back with gateway_mt and writes the ModEM,
WSINV3DMT, GoFEM and MARE2DEM files and the model with gateway_mt's writers,
so this catches any change in what those readers and writers produce.
"""
import importlib.util
import math
import os

from conftest import SYNTHETIC, quiet


def _load_generator():
    path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data', 'make_synthetic.py')
    spec = importlib.util.spec_from_file_location('make_synthetic', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _files(folder):
    found = []
    for root, _, names in os.walk(folder):
        found += [os.path.relpath(os.path.join(root, n), folder) for n in names]
    return sorted(f for f in found if not f.endswith('README.md'))


def _lines(path):
    with open(path, 'r') as f:              # text mode: CRLF / LF both read as '\n'
        return f.read().splitlines()


def _same_token(a, b):
    """Equal text, or numbers equal to 1 part in 1e9 (numpy builds can differ in the last digit)."""
    if a == b:
        return True
    try:
        return math.isclose(float(a), float(b), rel_tol=1e-9, abs_tol=1e-9)
    except ValueError:
        return False


def assert_same_file(new, old, name):
    new_lines, old_lines = _lines(new), _lines(old)
    assert len(new_lines) == len(old_lines), name
    for number, (a, b) in enumerate(zip(new_lines, old_lines), start=1):
        ta, tb = a.replace(',', ' ').split(), b.replace(',', ' ').split()
        assert len(ta) == len(tb) and all(map(_same_token, ta, tb)), \
            '{} line {}:\n  new: {}\n  old: {}'.format(name, number, a, b)


def test_regenerating_the_synthetic_data_reproduces_the_committed_files(tmp_path, monkeypatch):
    gen = _load_generator()
    out = tmp_path / 'synthetic'
    monkeypatch.setattr(gen, 'OUT', str(out))
    monkeypatch.setattr(gen, 'EDI_DIR', str(out / 'edis'))
    with quiet():
        gen.main()
    assert _files(out) == _files(SYNTHETIC)
    for name in _files(SYNTHETIC):
        assert_same_file(out / name, os.path.join(SYNTHETIC, name), name)
