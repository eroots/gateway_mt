"""Shared fixtures for the gateway_mt test suite.

Most tests run against the synthetic dataset in ``tests/data/synthetic``
(see its README). GUI tests use Qt's offscreen platform, so no display is
needed.
"""
import contextlib
import io
import os
import shutil

os.environ.setdefault('QT_QPA_PLATFORM', 'offscreen')
os.environ.setdefault('MPLBACKEND', 'Agg')

import pytest  # noqa: E402

import gateway_mt.data_structures as DS  # noqa: E402

SYNTHETIC = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data', 'synthetic')


def synthetic(*parts):
    """Path to a file in the synthetic dataset."""
    return os.path.join(SYNTHETIC, *parts)


@contextlib.contextmanager
def quiet():
    """Silence gateway_mt's progress printing (``Reading file: ...``)."""
    with contextlib.redirect_stdout(io.StringIO()):
        yield


@pytest.fixture(scope='session')
def dataset():
    """The full synthetic dataset: raw EDIs, inversion data, response and model.

    Session-scoped for speed, so tests must not modify it (use ``fresh_dataset``).
    """
    with quiet():
        return DS.Dataset(listfile=synthetic('edis', 'synthetic.lst'),
                          datafile=synthetic('synthetic_modem.dat'),
                          responsefile=synthetic('synthetic_modem_resp.dat'),
                          modelfile=synthetic('synthetic.model'))


@pytest.fixture
def fresh_dataset():
    """A private copy of the synthetic dataset that a test may modify."""
    with quiet():
        return DS.Dataset(listfile=synthetic('edis', 'synthetic.lst'),
                          datafile=synthetic('synthetic_modem.dat'),
                          responsefile=synthetic('synthetic_modem_resp.dat'))


@pytest.fixture
def edi_folder(tmp_path):
    """A temporary copy of the synthetic EDI folder (for commands that write next to the data)."""
    folder = tmp_path / 'edis'
    shutil.copytree(synthetic('edis'), folder)
    return folder


@pytest.fixture(scope='session')
def qapp():
    """A QApplication for GUI tests (offscreen)."""
    from PyQt5 import QtWidgets
    app = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])
    yield app
