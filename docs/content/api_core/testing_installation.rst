.. _Testing Installation:

Testing Installation
====================

A small synthetic dataset is included in the repository under
``tests/data/synthetic`` (30 sites; see the README in that folder for how it
was made). pip doesn't install it, so if you installed with pip, get a copy
of the repository first (or download it as a zip from GitHub)::

	git clone https://github.com/eroots/gateway_mt.git
	cd gateway_mt/tests/data/synthetic

From that folder, try launching each of the GUIs from the command line:

:ref:`Data Plot`::
 data_plot -n synthetic.pymt
:ref:`Mesh Designer`::
 mesh_designer synthetic_modem.dat
:ref:`Model Viewer`::
 model_viewer synthetic_modem.dat synthetic.model

If they all launch, you should be good to go.
If not, most likely there is a dependency issue somewhere.
:ref:`First Project` goes through a whole session with the same data.

From the top of the repository you can also run the automated tests, which
need no display::

	pip install -e .[test]
	pytest
