Getting Started
==========================

Launching the GUI
-----------------

The model viewer GUI is launched from the command line (or from the gateway_mt project window) via the command::

	model_viewer <model_file> <data_file>

The order of the files does not matter; they are recognised by their extensions (.model, .rho or .zani for models, and .dat or .data for data files).

It can take a while to load the GUI, particularly the first time you do so, due to the somewhat heavy dependencies required for 3-D viewing.

Required Files
--------------

The model plotting GUI requires at minimum a :ref:`Model File`, and an optional :ref:`Data File`, which is currently only used to plot station locations.

Default Behavior
----------------

Assuming valid files are used, the GUI should launch into the 2-D view, with the slice locations set to X=1, Y=1, Z=1 (bottom, left, and surface slices, respectively)

The 3-D view is initialized to a top-down (XY) view.

Default colour map is 'turbo' (inverted, so that conductors are red), with lower and upper colour map limits of 1 and 5 (log10 scale).

Currently the data file is used only to plot the station locations over the model. Plotting of induction arrows and phase tensors will likely be added in the future.


