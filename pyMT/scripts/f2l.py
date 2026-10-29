#!/usr/bin/env python
"""Build a ``.lst`` site-list file from the EDI / ``.dat`` files in the
current directory.

Usage::

    f2l              # list every .edi / .dat in the folder
    f2l simple       # same (explicit)
    f2l <suffix>     # only files whose name contains <suffix>

Installed as the ``f2l`` console script.
"""
import os
import sys

def search_dir(suffix):
    lst = []
    for file in os.listdir():
        print(file)
        if suffix in file and (file.endswith('edi') or file.endswith('.dat')):
            lst.append(file)
    return lst


def write_list(transect, lst, data_type):
    if data_type.lower() == 'amt':
        suffix = 'AMT'
        ending = 'A'
    elif data_type.lower() == 'bb':
        suffix = 'BB'
        ending = 'ML'
    elif data_type.lower() == 'all':
        suffix = 'all'
        ending = 'AML'
    else:
        suffix = ''
        ending = ''
    if data_type:
        lst = [file for file in lst if file[-5] in ending]
    else:
        lst = [file for file in lst if file.endswith(('edi', '.dat'))]
    lst = sorted(lst)
    with open('all'.join([transect, suffix, '.lst']), 'w') as f:

        f.write(str(len(lst)) + '\n')
        for file in lst:
            f.write(file + '\n')

def main():
    try:
        if sys.argv[1].lower() == 'simple':
            suffix = ['']
            data_types = ['']
        else:
            suffix = [sys.argv[1]]
            data_types = ['']
    except IndexError:
        suffix = ['']
        data_types = ['']
        # suffix = ['DRY', 'ATT', 'RRV', 'STU',
        #           'CHI', 'COB', 'GER', 'LAR',
        #           'MAL', 'MAT', 'ROU', 'SUD', 'SWZ', '']
        # data_types = ['amt', 'bb', 'all']
    for suff in suffix:
        lst = search_dir(suff)
        for data_type in data_types:
            write_list(suff, lst, data_type)


if __name__ == '__main__':
    main()    
