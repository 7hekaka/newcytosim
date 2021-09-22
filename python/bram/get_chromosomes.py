#!/usr/bin/env python3
#
# get_chromosomes.py
#
# simple extraction of data from files
#
# Copyright F. Nedelec, 2.11.2020

"""
    Collect data from given run directories, and print them to standard output

Syntax:

    get_chromosomes.py DIRECTORIES
    
Example:

    get_chromosomes.py run???? > data.txt

Description:

    This script is used to analyze chromosomes simulations

F. Nedelec, 2.11.2020
"""

import sys, os, subprocess, math

from pyned import find_differences, uncode, format_line
import read_config

import matplotlib.pyplot as plt

#------------------------------------------------------------------------


def find_indices(line):
    """
    return index after first 'nan'
    """
    s = len(line)
    a = 0
    while a < s and line[a] != 'nan':
        a += 1
    n = a+1
    while n < s and line[n] != 'nan':
        n += 1
    return a+1, n-1


def make_plot(arg):
    fts = 14
    fig = plt.figure(figsize=(8, 6))
    inx = find_indices(arg[0])[0]
    # plot each data point with a symbol:
    for data in arg:
        X = data[0]
        Y = data[inx+4]
        plt.plot(X, Y, 'o', markersize=4, markerfacecolor='none', markeredgecolor='blue')
    # label axes
    plt.xlabel('X', fontsize=fts)
    plt.ylabel('Y', fontsize=fts)
    #plt.xlim(3, 15)
    #plt.ylim(45, 360)
    plt.title('Simulations', fontsize=fts)
    fig.tight_layout()


#------------------------------------------------------------------------

def column_averages(file, start):
    """
    Calculate mean of values for each column, beyond line 'start'
    """
    inx = 0
    cnt = 0
    res = []
    from operator import add
    for line in file:
        vals = [ float(s) for s in uncode(line).split() ]
        inx = inx + 1
        if inx > start:
            cnt = cnt + 1
            if cnt > 1:
                res = list(map(add, res, vals))
            else:
                res = vals
    if cnt > 0:
        return [ round(v/cnt, 3) for v in res ]
    return [ 0 ]


#------------------------------------------------------------------------

def process(path):
    """ 
        This extracts parameters from the config file, and processes quantities
    """
    if path.startswith('run'):
        res = [ path[3:] ]
    else:
        res = [ path ]
    if not os.path.isfile(path+'/config.cym'):
        return []        
    try:
        par = find_differences('config.cym', path+'/config.cym')
    except IOError as e:
        par = 'unknown_parameters'
    res.extend(par)
    os.chdir(path)
    if not os.path.isfile('objects.cmo'):
        return []
    if not os.path.isfile('properties.cmo') and not os.path.isfile('properties.cmp'):
        return []
    res.append('nan')
    # get forces:
    filename = 'align.txt'
    if not os.path.isfile(filename):
        subprocess.call(['report2', 'chromosome:orientation', 'verbose=0'], stdout=open(filename, 'w'))
    with open(filename, 'r') as f:
        data = column_averages(f, 0)
        res.extend(data)
    return res

#------------------------------------------------------------------------

def main(args):
    paths = []
    
    for arg in args:
        if os.path.isdir(arg):
            paths.append(arg)
        else:
            sys.stderr.write("  Error: unexpected argument `%s'\n" % arg)
            sys.exit()
    if not os.path.isfile('config.cym'):
        sys.stderr.write("  Error: missing comparison base `config.cym'\n")
        sys.exit()
    
    if not paths:
        data = parse('.')
        print(format_line(data))
    else:
        res = []
        nb_columns = 0
        cdir = os.getcwd()
        for p in paths:
            #sys.stdout.write('- '*32+p+"\n")
            os.chdir(cdir)
            data = process(p)
            if not data:
                continue
            print(format_line(data))
            if nb_columns != len(data):
                if nb_columns == 0:
                    nb_columns = len(data)
                else:
                    sys.stderr.write("Error: data size mismatch in %s\n" % p)
                    break
            res.append(data)
        os.chdir(cdir)
        make_plot(res)
        plt.title('Simulations '+os.path.basename(os.getcwd()), fontsize=18)
        plt.savefig('simul.pdf', dpi=300)
        #plt.show()
        plt.close()


#------------------------------------------------------------------------

if __name__ == "__main__":
    if len(sys.argv) < 2 or sys.argv[1].endswith("help"):
        print(__doc__)
    else:
        main(sys.argv[1:])

