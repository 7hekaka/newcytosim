#!/usr/bin/env python3
#
# get_contraction.py
#
# simple extraction of data from files
#
# Copyright F. Nedelec, 2011 - 2016

"""
    Collect data from given run directories, and print them to standard output

Syntax:

    get_contraction.py DIRECTORIES
    
Example:

    get_contraction.py run???? > data.txt

Description:

    This script must be customized for any meaningful application,
    but it should be a useful template to start from.
    Please make a copy of the script with a different name.

F. Nedelec, Jan. 2016
"""

import sys, os, subprocess
from pyned import find_differences, uncode, get_column
import read_config

#------------------------------------------------------------------------

def get_parameters(path):
    pile = read_config.parse(path)
    res = {}
    for name in ['motor', 'complex', 'crosslinker']:
        try:
            cmd = read_config.get_command(pile, ['new', 'couple', name])
            res[name] = cmd.cnt
        except:
            res[name] = 0
    return res
    
    
#------------------------------------------------------------------------

def process(path):
    """ 
        This extracts parameters from the config file,
        and values from 'mom.txt'
    """
    if path.startswith('run'):
        res = path[3:] + ' '
    else:
        res = path + ' '
    try:
        if 0:
            nbs = get_parameters(path+'/config.cym')
            #res += '%8i %8i' % (2*nbs['motor']+nbs['complex'], 2*nbs['crosslinker']+nbs['complex'])
            res += '%8i %8i' % (nbs['motor'], nbs['crosslinker'])
        else:
            dif = find_differences('config.cym', path+'/config.cym')
            for x in dif:
                res += ' %9i' % x
    except IOError as e:
        sys.stderr.write("Error: %s\n" % repr(e))
        return ''
    res += ' nan'
    try:
        # get values from file 'mom.txt'
        with open(path+'/mom.txt', 'r') as f:
            res += get_column(f, -1)
    except IOError as e:
        sys.stderr.write("Error: %s\n" % repr(e))
        return ''
    return res


def main(args):
    paths = []
    
    for arg in args:
        if os.path.isdir(arg):
            paths.append(arg)
        else:
            sys.stderr.write("  Error: unexpected argument `%s'\n" % arg)
            sys.exit()
    
    if not paths:
        sys.stderr.write("  Error: you must specify directories\n")
        sys.exit()
    
    nb_columns = 0
    for p in paths:
        res = process(p)
        # check that the number of column has not changed:
        if nb_columns != len(res.split()):
            if nb_columns == 0:
                nb_columns = len(res.split())
            else:
                sys.stderr.write("Error: data size mismatch in %s\n" % p)
                return
        print(res)


#------------------------------------------------------------------------

if __name__ == "__main__":
    if len(sys.argv) < 2 or sys.argv[1].endswith("help"):
        print(__doc__)
    else:
        main(sys.argv[1:])


