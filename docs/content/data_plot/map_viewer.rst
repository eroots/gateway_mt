.. _Map Viewer:

Data Plot - Map Viewer
======================

.. figure:: ../../images/map_viewer.png
    :align: center
    :scale: 50 %

The Map Viewer is launched by clicking the :ref:`Show Map` button within the :ref:`Plot Options` tab.
This window is used to view the locations of the stations, as well as plot various data types in map view.
The window is broken into a `Menu Bar`_, `Control Dock`_, and the `Map Window`_.

Menu Bar
--------

.. figure:: ../../images/map_viewer_menu_bar.png
    :align: center
    :scale: 50 %

The Menu Bar has a number of options to customize the plotting within the `Map Window`_.

Colour Options
^^^^^^^^^^^^^^

The Colour Options menu contains options for colour map and colour limit selection.

Use the Colour Map menu to select the colour map to be used within the Map Window.

The Colour Limits menu is used to customize the lower and upper data limits to be plotted. Separate colour limits may be used for each type of data (e.g., apparent resistivity, phase, phase splits, skew, tipper, differences, and RMS).

Point Options
^^^^^^^^^^^^^

The Point Options menu has options for customizing the appearance of the point related data.

The Site Annotations sub-menu has options for controlling how stations are annotated.
By default, only the sites currently active within the :ref:`Plot Window` are annotated.
Annotations can also be turned off completely, or set to have all sites labelled.

The Site Markers sub-menu contains further menus for controlling the appearance of the site markers.

The Phase Tensors sub-menu is used the control size of the phase tensor ellipses, as well as the error tolerances (ellipses with larger errors are not plotted), the axis the azimuth is measured from, and the minimum axis ratio of the ellipses.

The Induction Arrows sub-menu is used to change the relative length of the induction arrows, the error tolerance and cutoff length, and the arrow convention (Parkinson (reversed), or Wiese).

The Coordinate System sub-menu sets whether stations are plotted in local, UTM, lat/long, or Lambert coordinates. Note that which of these is available will depend on what data is loaded (e.g., a ModEM data file alone has no information about the geographic locations of the stations).

The RMS Plotting sub-menu controls the size, colours and labels used by `Plot RMS`_.

Other Menus
^^^^^^^^^^^

The Axis Options menu sets the spatial units (km or m), an equal aspect ratio, and the font sizes. 'Lock Axis Limits' keeps the current zoom when the plot is redrawn.

The JPEG menu loads a geo-referenced JPEG image (with its world file) into the background of the map, and sets its opacity. The Coordinate System needs to be set to match the image.

The Pseudosection Options menu controls whether outliers are included in the `Pseudosection`_ interpolation.

Note that the scaling of the phase tensor ellipses and induction arrows has not be tested on all survey sizes, and so will likely need to be modified. The phase tensor ellipses are pre-normalized and so should generally plot well, however the induction arrows are plotted as is by default. In most cases, it is best to normalize their length (see the following section)

.. _Control Dock:

Control Dock
------------

.. figure:: ../../images/map_viewer_control_dock.png
    :align: center
    :scale: 50 %

The Control Dock is the main control panel for the Map Viewer. The various group boxes give options for plotting induction arrows, phase tensor ellipses, and apparent resistivity and phase pseudosection.

The currently plotted period / frequency is seen near the buttom of the Control dock, and can be changed using the nearby horizontal slider bar.

Induction Arrows
^^^^^^^^^^^^^^^^

The Induction Arrows group box allows for plotting of the induction arrows (in the convention chosen in the `Point Options`_ menu) within the `Map Window`_.

The Data and Response buttons plot the the induction arrows from the inversion data and response files in black and red, respectively. The drop-down menu selects real, imaginary, or both arrows.

The Normalize button scales the arrows so that they all have the same length. This is nearly always required in the current release, as a single noisy high amplitude induction arrow will tend to drown out all the others if the lengths are not normalized.

Phase Tensor
^^^^^^^^^^^^

The Phase Tensor groupbox is used to plot phase tensor ellipses in the `Map Window`_. The fill value of the ellipses is controlled by the contained drop-down menu.

The Data and Response checkboxes plot the phase tensor ellipses from the inversion data and response files, respectively.
If both checkboxes are selected, a phase tensor difference is plotted, using the definition chosen in the 'Difference Definition' drop-down menu. For example, the Heise definition (Heise et al.) is the misfit tensor

.. math ::
	\mathbf{\Delta} = \mathbf{I} - \frac{1}{2}(\mathbf{\Phi}_r^{-1}\mathbf{\Phi}_d + \mathbf{\Phi}_d\mathbf{\Phi}_r^{-1})

where :math:`\mathbf{\Phi}_d` and :math:`\mathbf{\Phi}_r` are the data and response phase tensors, while the Booker definition is simply :math:`\mathbf{\Phi}_d - \mathbf{\Phi}_r`.
In this case, the colour of the ellispes may also be filled by the :math:`\delta` value give as a percentage, defined by:

.. math ::
	\delta = 100 * \left\| \mathbf{\Phi}_d - \mathbf{\Phi}_r \right\| / \left\| \mathbf{\Phi}_d \right\|

and the other fill values give the difference in that parameter (response minus data).

All other phase tensor fill values are defined as in Caldwell et al., 2004. This includes the skew value :math:`\beta`, and so the recommended upper limit for approximate two-dimensionality is :math:`\beta \leq \left|3^{\circ} \right|`

Pseudosection
^^^^^^^^^^^^^

The Pseudosection groupbox gives options for plotting map view sections of the apparent resistivity and phase, as well as other data (Bostick depth, tipper amplitude, and phase tensor parameters). The sections are interpolated with a radial basis function, linear, cubic, or nearest-neighbour scheme. If the naturalneighbor python package (https://pypi.org/project/naturalneighbor/) is installed, the Natural Neighbor interpolation scheme of Sibson (1981) is also available.

The first drop-down menu controls what is plotted.

The Data and Response checkboxes plot the sections from the inversion data or response files, respectively, while checking both boxes will result in difference plots. The difference in apparent resistivities is expressed as a percent difference, while for phase it is given in degrees.

The XY, YX, and determinant resistivities and phases (and the XY-YX difference) may be toggled through the second drop-down menu.

The Interp points spinbox controls how many grid points are used in the section interpolation. The default is 200 points. Few points may be used to increase the responsiveness of the GUI, while more points may be used for a finer grid.

Plan View
^^^^^^^^^

If a :ref:`Model File` is loaded, the Plan View groupbox plots a depth slice of the model underneath the stations.

.. _Plot RMS:

Plot RMS
^^^^^^^^

The relative RMS misfit of each station can be overlaid by checking the Plot RMS checkbox. Unlike the other data plotting options in this window, the RMS misfit overlay does not change per period, but reflects instead the total RMS misfit across all periods and components for each site.

.. _Map Window:

Map Window
----------

.. figure:: ../../images/map_viewer_map_window.png
    :align: center
    :scale: 50 %

The Map Window is the plot area of the Map Viewer. The site locations are plotted here, along with any data components that have been activated in the `Control Dock`_.

By default, only the station locations are plotted along with the annotation style set in the `Point Options`_ menu. A colour bar will be added if necessary, e.g., when plotting phase tensor ellipses.

Note that by default the subplot used is set to fill the available space, and therefore does not have an equal aspect ratio (see `Other Menus`_).

The toolbar at the bottom of the Map Window may be used to pan, zoom, and customize the plot itself.

The Map Window is reset any time an element of the plot changes. As such, zoom and pans that are applied will not be remembered when the plot is changed, unless 'Lock Axis Limits' is checked (see `Other Menus`_).
