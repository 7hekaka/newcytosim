#!/usr/bin/env python
#
# get_polymer.py
#
# simple extraction of data from files
#
# Copyright F. Nedelec, 20.02.2017

"""
    Collect data from given run directories, and print them to standard output

Syntax:

    get_polymer.py DIRECTORIES
    
Example:

    get_polymer.py run???? > data.txt

Description:

    This script was customized for Platelet Simulations.

F. Nedelec, Feb. 2017
"""

import sys, os, subprocess
import read_config
from pyned import find_differences, uncode

#------------------------------------------------------------------------

def get_parameters(path):
    res = []
    pile = read_config.parse(path)
    cmd = read_config.get_command(pile, ['set', 'fiber', 'microtubule']);
    pam = cmd.values()
    try:
        val = pam['hydrolysis_rate']
    except:
        val = find_differences('config.cym', path+'/config.cym')
    return [val]


def get_values(path):
    cdir = os.getcwd();
    os.chdir(path)
    sub = subprocess.Popen(['report3', 'fiber:length'], stdout=subprocess.PIPE)
    # Get results from standard output:
    res = 0
    cnt = 0
    mts = 0
    for stuff in sub.stdout:
        line = uncode(stuff).split()
        #print(line)
        if len(line) < 2:
            pass
        elif line[0] == '%':
            if line[1] == 'start':
                sec = float(line[2])
            elif line[1] == 'end':
                pass
        elif len(line) == 7:
            mts += float(line[1])
            res += float(line[6])
            cnt += 1
    sub.stdout.close()
    os.chdir(cdir);
    return (mts/cnt, res/cnt)


#------------------------------------------------------------------------

def process(path):
    """ 
        This extracts parameters from the config file,
        and values from 'mom.txt'
    """
    # put file name:
    if path.startswith('run'):
        res = path[3:] + ' '
    else:
        res = path + ' '
    # add parameters:
    pam = get_parameters(path+'/config.cym')
    for v in pam:
        res += ' ' + repr(v)
    res += ' nan'
    # add data obtained with 'report'
    vals = get_values(path)
    for v in vals:
        res += ' %12.3f' % v
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


