#!/usr/bin/env python
"""Build a ``.lst`` site-list file from the EDI / j-format files in a folder.

Usage::

    f2l                 # every .edi / .dat site -> all_sites.lst
    f2l simple          # same as no filter
    f2l swz             # only sites whose file name contains 'swz' -> swz_sites.lst
    f2l -d edis -o my_sites.lst --force

Site names are written without their extension (``read_sites`` does not need
them). The filter is case-insensitive. Installed as the ``f2l`` console script.
"""
import argparse
import os
import sys

EXTENSIONS = ('.edi', '.dat')


def search_dir(path, name_filter=''):
    """Site names (no extension) of the .edi / .dat files in ``path`` matching ``name_filter``."""
    sites = []
    for file in os.listdir(path):
        name, ext = os.path.splitext(file)
        if ext in EXTENSIONS and name_filter.lower() in file.lower():
            if name not in sites:   # SITE.edi + SITE.dat -> one entry (gateway_mt reads the .dat)
                sites.append(name)
    return sorted(sites)


def write_list(sites, outfile):
    with open(outfile, 'w') as f:
        f.write(str(len(sites)) + '\n')
        for site in sites:
            f.write(site + '\n')


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        prog='f2l',
        description='Write a gateway_mt site list (.lst) of the EDI / j-format (.dat) files in a folder.',
        epilog='Note: every .dat file counts as a site, so keep inversion data files '
               '(e.g. ModEM .dat) out of the folder or use a FILTER.')
    parser.add_argument('filter', nargs='?', default='',
                        help="only include files whose name contains this text (case-insensitive); "
                             "omit it or use 'simple' to include every site")
    parser.add_argument('-d', '--dir', default='.',
                        help='folder to search (default: current folder)')
    parser.add_argument('-o', '--output',
                        help="output file (default: all_sites.lst, or <FILTER>_sites.lst with a filter); "
                             "written inside --dir unless a path is given")
    parser.add_argument('-f', '--force', action='store_true',
                        help='overwrite the output file if it already exists')
    args = parser.parse_args(argv)
    if args.filter.lower() == 'simple':
        args.filter = ''
    return args


def main(argv=None):
    args = parse_args(argv)
    if not os.path.isdir(args.dir):
        print('f2l: {} is not a folder'.format(args.dir))
        return 1
    sites = search_dir(args.dir, args.filter)
    if not sites:
        print('f2l: no .edi / .dat files{} in {}'.format(
            " matching '{}'".format(args.filter) if args.filter else '', os.path.abspath(args.dir)))
        return 1
    outfile = args.output or '{}_sites.lst'.format(args.filter or 'all')
    if not os.path.dirname(outfile):
        outfile = os.path.join(args.dir, outfile)
    if os.path.exists(outfile) and not args.force:
        print('f2l: {} already exists (use --force to overwrite)'.format(outfile))
        return 1
    write_list(sites, outfile)
    print('Wrote {} site{} to {}'.format(len(sites), '' if len(sites) == 1 else 's', outfile))
    return 0


if __name__ == '__main__':
    sys.exit(main())
