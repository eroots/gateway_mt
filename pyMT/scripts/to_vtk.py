"""Write a model and / or station locations to VTK files for ParaView.

Two ways to run it:

* ``to_vtk`` with no arguments asks for the inputs interactively (model
  file, list or data file, coordinate system, output name).
* ``to_vtk -o NAME [options]`` does the same from the command line; see
  ``to_vtk -h``. For example::

      to_vtk -o survey -m inv.model -l all_sites.lst --utm 16U

Writes ``NAME_model.vtk`` and / or ``NAME_sites.vtk``. Meant for ModEM /
WSINV3DMT-style models -- GoFEM already writes its models as VTU.
Installed as the ``to_vtk`` console script. See
``examples/to_vtk_point-cloud.py`` for a point-cloud variant.
"""
import argparse
import os
import sys

import pyproj

import pyMT.data_structures as WSDS
from pyMT.IO import verify_input
from pyMT.WSExceptions import WSFileError
from pyMT.utils import project


def transform_locations(dataset, UTM):
    """Re-project the raw-data site locations to the given UTM zone (or Lambert)."""
    dataset.raw_data.locations = dataset.raw_data.get_locs(mode='latlong')
    if 'lam' in UTM.lower():
        print("Reminder: Current Lambert Transformation is set to EPSG3978")
        transformer = pyproj.Transformer.from_crs('epsg:4326', 'epsg:3978')
        for ii, (lat, lon) in enumerate(dataset.raw_data.locations):
            x, y = transformer.transform(lat, lon)
            dataset.raw_data.locations[ii, :] = y, x

    else:
        UTM_letter = UTM[-1]
        while True:
            if len(UTM) == 3:
                UTM_number = int(UTM[:2])
                break
            elif len(UTM) == 2:
                UTM_number = int(UTM[0])
                break
            else:
                print('UTM {} is not a valid zone.'.format(UTM))
                UTM = verify_input('UTM Zone',
                                   expected=str,
                                   default='dummy')
        for ii in range(len(dataset.raw_data.locations)):
            easting, northing = project((dataset.raw_data.locations[ii, 1],
                                         dataset.raw_data.locations[ii, 0]),
                                        zone=UTM_number, letter=UTM_letter)[2:]
            dataset.raw_data.locations[ii, 1], dataset.raw_data.locations[ii, 0] = easting, northing


def to_vtk(outfile, datafile=None, listfile=None, modelfile=None,
           datpath=None, origin=None, UTM=None, sea_level=0, use_elevation=True,
           resolutionfile=None, transform_coords=None, trim=None):
    """Write ``outfile_model.vtk`` and / or ``outfile_sites.vtk``.

    Site locations come from ``listfile`` (raw data) if given, else from
    ``datafile``. When ``transform_coords == 'y'`` the raw-data sites are
    re-projected to ``UTM`` and the origin is recomputed from them, so any
    ``origin`` passed in is ignored.
    """
    if not outfile:
        print('Output file required!')
        return
    if listfile or datafile:
        dataset = WSDS.Dataset(modelfile=modelfile,
                               datafile=datafile,
                               listfile=listfile,
                               datpath=datpath)
    if not UTM:
        print('Using dummy UTM zone')
        UTM = '999'
    if transform_coords == 'y':
        if listfile:
            print('Transforming locations to UTM {}'.format(UTM))
            transform_locations(dataset, UTM)
        else:
            print('Coordinate transformation needs a list file (raw data); skipping it.')
            transform_coords = 'n'

    if not origin and not listfile:
        print('You must either specify the origin, or the list file to calculate it from.')
        return
    elif (not origin and listfile) or transform_coords == 'y':
        if origin and transform_coords == 'y':
            print('Ignoring the given origin; using the centre of the transformed sites.')
        origin = dataset.raw_data.origin
        print('Setting origin as ({:.1f}, {:.1f})'.format(*origin))
    if modelfile:
        model = WSDS.Model(modelfile=modelfile)
        model.origin = origin
        model.UTM_zone = UTM
        print('Writing model to {}'.format('_'.join([outfile, 'model.vtk'])))
        if trim:
            for ix in range(trim[0]):
                model.dx_delete(0)
                model.dx_delete(model.nx)
            for ix in range(trim[1]):
                model.dy_delete(0)
                model.dy_delete(model.ny)
            for ix in range(trim[2]):
                model.dz_delete(model.nz)
        if resolutionfile:
            print('Adding resolution')
            resolution = WSDS.Model(modelfile=resolutionfile)
            model.resolution = resolution.vals
        model.to_vtk(outfile, sea_level=sea_level)
    if listfile:
        print('Writing sites to {}'.format('_'.join([outfile, 'sites.vtk'])))
        dataset.raw_data.locations -= (origin[1], origin[0])
        dataset.raw_data.to_vtk(origin=origin, UTM=UTM, outfile=outfile,
                                sea_level=sea_level, use_elevation=use_elevation)
    elif datafile:
        print('Writing sites to {}'.format('_'.join([outfile, 'sites.vtk'])))
        dataset.data.to_vtk(origin=origin, UTM=UTM, outfile=outfile,
                            sea_level=sea_level, use_elevation=use_elevation)


def is_gofem(path):
    if path and path.lower().endswith('.gdat'):
        print('{} is GoFEM data -- to_vtk is for ModEM / WSINV3DMT-style inputs '
              '(GoFEM writes its own VTU files).'.format(path))
        return True
    return False


def get_inputs():
    """Ask for the inputs interactively; returns keyword arguments for :func:`to_vtk`."""
    raw_origin = '0, 0'

    args = {}
    to_output = verify_input(message='Output model, sites, or both? (m/d/b)',
                             expected='mbd', default='b')
    if not to_output:
        to_output = 'b'
    if to_output == 'm' or to_output == 'b':
        # Get model name
        modelfile = verify_input('Input model name', expected='read')
        args.update({'modelfile': modelfile})
        include_resolution = verify_input('Include model resolution?',
                                          expected='yn', default='n')
        if include_resolution == 'y':
            resolution_file = verify_input('Resolution file name:', expected='read',
                                           default='Resolution0_inverted.model')
            args.update({'resolutionfile': resolution_file})
        trim = verify_input('Trim model?', expected='yn', default='n')
        if trim == 'y':
            trim_x = verify_input('Input trim (north-south)', expected=int, default=5)
            trim_y = verify_input('Input trim (east-west)', expected=int, default=5)
            trim_z = verify_input('Input trim (depth)', expected=int, default=5)
            args.update({'trim': (trim_x, trim_y, trim_z)})
    if to_output == 'd' or to_output == 'b':
        # Get data or list file
        while True:
            datafile = verify_input('Input data or list file name', expected='read')
            if not is_gofem(datafile):
                break
        if datafile.lower().endswith('.lst'):
            # If listfile, make sure the conversion to raw data works
            path = os.path.dirname(os.path.abspath(datafile))
            datpath = input('Path to EDI / .dat files \n' +
                            '[Default: {}] > '.format(path)) or path
            try:
                raw_data = WSDS.RawData(listfile=datafile, datpath=datpath)
                raw_origin = '{}, {}'.format(*[float(x) for x in raw_data.origin])
                args.update({'datpath': datpath})
                args.update({'listfile': datafile})
            except WSFileError as err:
                print(err.message)
        else:
            args.update({'datafile': datafile})
    args.update({'UTM': verify_input('UTM Zone',
                                     expected=str,
                                     default='dummy')})
    if args['UTM'] != 'dummy' and 'listfile' in args:
        args.update({'transform_coords': verify_input('Do you want to transform the coordinates to this UTM zone?',
                                                      expected='yn',
                                                      default='y')})
    else:
        args.update({'transform_coords': 'n'})
    if args['transform_coords'] != 'y':
        # (When transforming, the origin is taken from the transformed sites.)
        origin = verify_input('Input model origin as Easting, Northing',
                              expected='numtuple',
                              default=raw_origin)
        args.update({'origin': origin})
    sea_level = verify_input('Sea level adjustment (Positive for above sea level)',
                             expected=float, default=0)
    args.update({'sea_level': sea_level})
    use_elevation = verify_input('Use site elevation?',
                                 expected='yn', default='y')
    args.update({'use_elevation': use_elevation == 'y'})
    outfile = verify_input('Base output file name', expected='write')
    args.update({'outfile': outfile})

    return args


def parse_args(argv):
    parser = argparse.ArgumentParser(
        prog='to_vtk',
        description='Write a model and / or station locations to VTK for ParaView. '
                    'Run with no arguments to be prompted for each input instead.',
        epilog='For ModEM / WSINV3DMT-style models; GoFEM already writes VTU.')
    parser.add_argument('-o', '--output', required=True,
                        help='base output name: writes OUTPUT_model.vtk and / or OUTPUT_sites.vtk')
    parser.add_argument('-m', '--model', help='model file')
    sites = parser.add_mutually_exclusive_group()
    sites.add_argument('-l', '--list', help='site list (.lst) -- sites from the raw EDI / .dat files')
    sites.add_argument('-d', '--data', help='inversion data file -- sites from the data file')
    parser.add_argument('--datpath', help='folder holding the EDI / .dat files (default: the list file\'s folder)')
    parser.add_argument('--utm', help="UTM zone, e.g. 16U (or 'lambert'); with --list, sites are "
                                      "re-projected to it unless --no-transform")
    parser.add_argument('--no-transform', action='store_true',
                        help='do not re-project the sites to --utm')
    parser.add_argument('--origin', nargs=2, type=float, metavar=('EASTING', 'NORTHING'),
                        help='model origin (default: centre of the sites from --list, else 0 0)')
    parser.add_argument('--sea-level', type=float, default=0,
                        help='sea-level adjustment, positive above sea level (default 0)')
    parser.add_argument('--no-elevation', action='store_true',
                        help='place all sites at the sea-level height instead of their elevations')
    parser.add_argument('--resolution', help='resolution file to add to the model')
    parser.add_argument('--trim', nargs=3, type=int, metavar=('NX', 'NY', 'NZ'),
                        help='cells to trim from each side in x / y and from the bottom in z')
    args = parser.parse_args(argv)
    if not (args.model or args.list or args.data):
        parser.error('give a model (-m) and / or sites (-l or -d)')
    return args


def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]
    if not argv:
        try:
            args = get_inputs()
        except KeyboardInterrupt:
            print('\nQuitting...')
            return 1
        to_vtk(**args)
        return 0
    args = parse_args(argv)
    if is_gofem(args.data):
        return 1
    transform = 'y' if (args.utm and args.list and not args.no_transform) else 'n'
    if args.origin is None and not args.list:
        print('No --origin given; using 0, 0.')
        args.origin = (0.0, 0.0)
    to_vtk(outfile=args.output, modelfile=args.model, listfile=args.list, datafile=args.data,
           datpath=args.datpath or (os.path.dirname(os.path.abspath(args.list)) if args.list else None),
           origin=tuple(args.origin) if args.origin else None, UTM=args.utm,
           sea_level=args.sea_level, use_elevation=not args.no_elevation,
           resolutionfile=args.resolution, transform_coords=transform,
           trim=tuple(args.trim) if args.trim else None)
    return 0


if __name__ == '__main__':
    sys.exit(main())
