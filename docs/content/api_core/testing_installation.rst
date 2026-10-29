.. _Testing Installation:

Testing Installation
====================

A small synthetic dataset is included in the repository under
``tests/data/synthetic`` (30 sites; see the README in that folder for how it
was made). From that folder, try launching each of the GUIs from the
command line:

:ref:`Data Plot`::
 data_plot -n synthetic.pymt
:ref:`Mesh Designer`::
 mesh_designer synthetic_modem.dat
:ref:`Model Viewer`::
 model_viewer synthetic_modem.dat synthetic.model

If they all launch, you should be good to go.
If not, most likely there is a dependency issue somewhere.
