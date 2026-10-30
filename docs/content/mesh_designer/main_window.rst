.. _Mesh Designer:

Mesh Designer - Main Window
===========================

.. figure:: ../../images/mesh_designer.png
    :align: center
    :scale: 50 %

The Mesh Designer has 3 areas: The `Menu Bar`_, `Control Dock`_, and `Plot Window`_.

Menu Bar
--------

.. figure:: ../../images/mesh_designer_menu_bar.png
    :align: center
    :scale: 50 %

The Save / Revert menu has options for writing the model to a file, as well as saving and reverting progess within the Mesh Designer.

The Write Model sub-menu will open a dialog box that asks for the name of the file to be written to. You will be prompted again if this operation is going to overwrite an existing file. A second prompt will ask for a name to output an associated covariance file (hit cancel to skip).

* The model can be written as a ModEM 3-D file, or as an MT3DANI (anisotropic) file. The ModEM format is the same as the WSINV3DMT file format, except that the resistivity values will be given as the natural logarithm.

* If you need the file in WSINV3DMT format, the workaround at present is to change the file format programmatically (see :ref:`Recipes`)

The Save Progress sub-menu internally saves any changes you've made to the mesh. This checkpoint can then be reverted to later by using the Revert Progress sub-menu.

The Colours menu is controls the colour map and colour limits used, as well as the colour of the mesh lines.

* As this GUI is meant mainly for editing the mesh, these options are mainly to ensure good visibility of the mesh lines regardless of the resistivity used

The Display Options menu controls the aspect ratio of the plot, and the appearance and annotation of the station markers.

The 1D Modelling menu opens the 1-D modelling window (as in :ref:`Data Plot`), whose synthetic or inverted model can be used as the starting model (see below).

Control Dock
------------

.. figure:: ../../images/mesh_designer_control_dock.png
    :align: center
    :scale: 50 %

The Control Dock is the main control panel for manipulating the mesh and model. It is broken into 3 tabs:

* `Manipulate Mesh`_
* `Reset Mesh`_
* `Smoothing`_

Manipulate Mesh
^^^^^^^^^^^^^^^

The Manipulate Mesh tab is used to, as the name suggests, to modify and manipulate the mesh used.

Specifically, there are 3 things that can be modified from this tab: The XY Padding, the Z mesh (or depth mesh), and the starting model.

The XY Padding groupbox is used to add and remove padding cells in the XY plane.

The Add Pad and Remove Pad buttons add and remove pads from the left / right / top bottom/ of the XY plane, depending on which of the corresponding checkboxes are selected.

When adding pads, the size of the new pad is determined by taking the size of the outer most cell and multiplying it by the value in the Pad Size Multiplier spinbox.

* Note: As the definition of the mesh within the ModEM and WSINV3DMT file formats has no explicit origin, all models generated here will have their origins placed in the center of the mesh.
* For this reason, it is important to ensure that the model is even on the left / right and top / bottom.
	* If it is not, the definition of the mesh relative to the station locations may not be the same as appears while using this GUI.
* Always double check the output mesh and data files with another tool afterwards (e.g., with :ref:`Model Viewer`)

The depth mesh is controlled through the Z Mesh groupbox, which has two ways of generating it: By Decade, and By Increase Factor.

In the By Decade tab, the thickness of the first slice (in meters) is specified in the Min. Depth edit line. The final depth (i.e., the maximum depth to use in the mode) is specified in the Max. Depth edit line.

Specification of the mesh between the first and last depth is controlled in the Depths / Decade list. This list will be automatically populated with a list of values. The length of this list is such that there is one value per decade of depth.

Once each of these values has been specified, hit the Generate Depths button to generate the Z mesh.

For example, for a first depth of 1 m and a last depth of 500000 m, the Depths / Decade list will be populated with 6 values. From top to bottom, they correspond to the number of layers used between depths of:

* 1-10 m
* 10-100 m
* 100-1000 m
* 1000-10000 m
* 10000-100000 m
* 100000-500000 m

In this instance, each value corresponds to the number of logarithmically spaced layers to use within each decade.

In the By Increase Factor tab, the mesh is instead specified by the factor by which each layer is thicker than the one above it. Each row of the Depth Cutoff / Increase Factors table gives a depth, and the factor to use down to that depth.

As a general rule of thumb, it is best to ensure that the sizes of the layers are always increasing. In accordance with this, the program will check the 2nd derivative of the generated mesh. If the derivative is negative anywhere, a message will appear in the Messages box saying so.

The Starting Model controls set the resistivity of the model. The backround resistivity of the model may be changed by editing the corresponding line and clicking the Set Background button. Alternatively, if the 1-D modelling window is open, the synthetic or inverted 1-D model can be used ('Use 1D Synthetic' / 'Use 1D Inverse') instead of a halfspace.

Reset Mesh
^^^^^^^^^^

The Reset Mesh tab is used to generate a new, uniformly spaced mesh from scratch.

Set the minimum cell sizes for the X and Y directions in the corresponding boxes, and hit Regenerate Mesh.

A new, uniform mesh will be generated with cells of the minimum size (the maximum sizes are not used here), extending to the bounds set by the station locations. The boxes start at 1/60 of the spread of the stations (at least 1 km). If the new mesh would be very large (more than 10 million cells), you will be asked before it is built.

Smoothing
^^^^^^^^^

The Smoothing tab is used to smooth the resistivity values of an inverted model.

* This tab is somewhat experimental at the moment *

Set the smoothing length in the X, Y, and Z direction using the Sigma X / Y / Z spinboxes, and hit Smooth Model to apply a Gaussian smoother with those parameters. The Z value in the Slices groupbox chooses which depth slice is shown in the `Plot Window`_.

Plot Window
-----------

.. figure:: ../../images/mesh_designer_plot_window.png
    :align: center
    :scale: 50 %

The Plot Window shows the current mesh, and if a data file was included, the station locations.

Fine grained modification of the mesh is done within this window.

Left click anywhere within the mesh to add a new vertical mesh line. Right click to add a horizontal mesh line.

Holding CTRL while left or right clicking will remove the nearest vertical / horizontal mesh line, respectively.

Refrain from double clicking within this window. Although some precautions have been implemented to avoid generated invalid meshes, double clicking can sometimes result in multiple mesh lines in the same location, i.e., a cell with 0 width, which will subsequently crash ModEM.

At the bottom of the Plot window is a toolbar. From here, you can pan and zoom into the plot, as well as return to the home view using the corresponding buttons.

Note that panning and zooming *is* stored within this window, which allows you to zoom into an area of high site density and add additional mesh:

* Click the zoom button and draw a rectangle around the area of interest.
* Unclick the zoom button to return the click functionality back to mesh modification.
* Modifiy the mesh as required.
* Hit the Home button in the toolbar to return to your original (un-zoomed) view.