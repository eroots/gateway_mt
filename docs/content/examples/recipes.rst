.. _Recipes:

Recipes
=======

All the following recipes assume gateway_mt modules have been imported as follows::

    import gateway_mt.data_structures as ds
    import gateway_mt.utils as utils

Adding Data to an Existing Data File
------------------------------------

Adding new periods to an existing data file is supported through the :ref:`Data Plot` GUI. However, it is assumed that the list file and data file used are consistent, i.e., that they have the same stations and in the same order.
Therefore, adding new stations to an existing data file must be done using the gateway_mt API.

This example uses the :ref:`Data`, and :ref:`RawData` classes as well as the RawData.get_data and Data.add_site methods.

Assuming the data file we want to modify is 'demo.dat', and the list file we want to pull in new sites from is 'raw/all.lst'::

    data = ds.Data('demo.dat')  # Load in the data file
    raw_data = ds.RawData('raw/all.lst')  # Load in the raw data
    # Get the new site to add with the same periods and components as 'data'
    to_add = raw_data.get_data(periods=data.periods,
                               components=data.components,
                               sites='test_1')
    data.add_site(to_add['test_1'])  # Add the new site to the data object
    data.write('new_data', file_format='modem')  # Write the new data out.

Alternatively, multiple ModEM data files can be combined (e.g., append tipper data from one file into another), so long as resulting file still adheres to gateway_mt's requirements (i.e., uniform periods and data components for all sites).

Working with Inconsistent Periods / Components
----------------------------------------------

Sometimes it may be desireable to invert a data set that has a non-uniform period set and / or component set, for example when mixing AMT and BB stations. If the number of wasted data points is low (maybe 3-4 frequencies total), it is usually simpler to just set high errors to the unavailable frequencies and let the inversion run as normal (the number of processors required is the same regardless, addition RAM required for wasted data points is low).
Note, a 'wasted data point' is a data point which is will be inverted at a subset of stations, but is unavailable at other stations (most often when mixing stations with different bandwidths, or stations with / without tipper)
However, if you intend to invert many frequencies and / or components that would be unavailable at a subset of stations, there is a workaround. This workaround requires some knowledge of Python.

First, you should set up the ModEM data file as usual (i.e., using Data Plot).
Make note of (or intentionally set up) any patterns you can use to distinguish your stations. In this example, broadband stations end in 'M', and AMT stations end in 'A'. The final 4 periods of AMT stations, and first 3 periods of BB stations will be flagged and removed::

    data = ds.Data('demo.dat')  # Load your data set
    for site in data.site_names:  # Loop through stations
        if site.lower().endswith('a'):  # If its an AMT station
            for comp in data.components:  # Loop over available components
                data.sites[site].used_error[comp][-4:] = data.REMOVE_FLAG  # Flag the errors for the last 4 periods
        elif site.lower().endswith('m'):
            for comp in data.components:
                data.sites[site].used_error[comp][:3] = data.REMOVE_FLAG  # Flag the errors for the first 3 periods
    data.write('demo_flagged.dat', write_removed=True)  # Write out the flagged and removed versions

This writes two files. 'demo_flagged.dat' has a uniform period band, with the flagged data written at the end of the file with errors equal to Data.REMOVE_FLAG (at the time of writing, it is 1234567). 'demo_flagged_removed.dat' is the same file with the flagged data left out.

You can now invert 'demo_flagged_removed.dat'. gateway_mt can read this file (and the response file created by inverting it) directly: any data missing from a station is filled in as dummy data and flagged, so that the period band is uniform again.
