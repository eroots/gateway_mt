Scripts
=======

A number of scripts come with pyMT, although only a few of these are ready for `General Use`_ . Others will need to be manually edited for your particular needs. Only a few will be described here.

General Use
-----------

The scripts listed here are ready to run out of the box, and either have command line inputs to control them, or are so simple as to not require any inputs.

to_vtk
^^^^^^

Script for converting site locations and models to VTK files (compatible with, e.g., Paraview).

Usage::

	to_vtk.py

Follow the command line prompts to select your model and/or data files. You will also have the option to project your model and data into a particular UTM zone.

Note that any station / model rotations will have to be handled manually.

autogen_report
^^^^^^^^^^^^^^

Writes a misfit report for an inversion result: a text summary (input files, inversion type, overall and per-component RMS, worst-fitting sites, the inverted periods, and blank fields for your own notes) and a PDF of data vs. response curves, four sites per page. Each site shows apparent resistivity and phase, or the impedance tensor with ``--plot impedance``, with the tipper underneath when present.

Usage::

	autogen_report -d <data file> -r <response file> [-l <list file>] [-m <model file>] [-o <report folder>] [--plot rho-phase|impedance]

Any data format pyMT can read is accepted (ModEM, WSINV3DMT, GoFEM, ...). Give the list file to add the raw data to the plots, and to name the sites when the data format doesn't store names (WSINV3DMT).

plot_dataset
^^^^^^^^^^^^

Plots raw data, inversion data and / or an inversion response -- any combination -- to a PDF, four sites per page, in the same layout as autogen_report. Useful for checking raw data, or an inversion data file before inverting.

Usage::

	plot_dataset [-l <list file>] [-d <data file>] [-r <response file>] [-o <output.pdf>] [--plot rho-phase|impedance] [--sites SITE ...]

