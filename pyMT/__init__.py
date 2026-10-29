"""Compatibility alias: ``pyMT`` was renamed ``gateway_mt`` in version 0.12.0.

``import pyMT.data_structures`` and friends keep working (they load the
``gateway_mt`` modules, so ``pyMT.utils is gateway_mt.utils``), but warn.
Change imports to ``gateway_mt``; this alias will be removed in a future
release.
"""
import importlib
import importlib.abc
import importlib.util
import sys
import warnings

warnings.warn("The pyMT package has been renamed gateway_mt. Please change 'import pyMT...' to "
              "'import gateway_mt...'; the pyMT name will be removed in a future release.",
              FutureWarning, stacklevel=2)


class _AliasFinder(importlib.abc.MetaPathFinder, importlib.abc.Loader):
    """Resolve ``pyMT.<name>`` to the already-renamed ``gateway_mt.<name>`` module."""
    prefix = 'pyMT.'

    def find_spec(self, fullname, path, target=None):
        if fullname.startswith(self.prefix):
            return importlib.util.spec_from_loader(fullname, self)
        return None

    def create_module(self, spec):
        return importlib.import_module('gateway_mt.' + spec.name[len(self.prefix):])

    def exec_module(self, module):
        pass


sys.meta_path.insert(0, _AliasFinder())

from gateway_mt import __version__  # noqa: E402,F401
