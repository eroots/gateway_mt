Installation
============

gateway_mt needs Python 3.10 or newer. The Python dependencies (numpy, scipy,
matplotlib, PyQt5, vtk, pyvista, pyvistaqt, pyproj, pyshp, colorcet, cma) are
installed automatically.

Installing gateway_mt
---------------------

The tested route uses Anaconda (or Miniconda):

1. Create and activate an environment::

	conda create -n gateway_mt python=3.12 pip
	conda activate gateway_mt

2. Install gateway_mt::

	pip install git+https://github.com/eroots/gateway_mt.git

3. To update later::

	pip install --upgrade git+https://github.com/eroots/gateway_mt.git

If you intend to modify the code, clone the repository and install it in
editable mode instead::

	git clone https://github.com/eroots/gateway_mt.git
	cd gateway_mt
	pip install -e .

Then check the installation with the bundled test data (see
:ref:`Testing Installation`).

Optional: natural-neighbour interpolation
-----------------------------------------

The plan-view pseudosections in :ref:`Data Plot` use scipy-based
interpolation by default. If the ``naturalneighbor`` package is installed, a
natural-neighbour scheme is also offered. An updated version that works with
current Python and numpy can be installed with::

	pip install git+https://github.com/eroots/natural-neighbor-interpolation

It needs the Microsoft Visual C++ build tools on Windows.

Coming from pyMT
----------------

gateway_mt is the continuation of pyMT under a new name; pyMT itself is still
available. The two are separate packages and install the same command names
(``gateway_mt``, ``data_plot``, ...), so install gateway_mt in its own
environment rather than alongside pyMT. In your own scripts, replace
``import pyMT...`` with ``import gateway_mt...``.
