#!/usr/bin/env python3
#
# get_ring.py
#
# simple extraction of data from files
#
# Copyright F. Nedelec, 2011 - 2016

"""
    Collect data from given run directories, and print them to standard output

Syntax:

    get_ring.py DIRECTORIES
    
Example:

    get_ring.py run???? > data.txt

Description:

    This script is used to analyze platelet simulations

F. Nedelec, July 2017, 25--27 Aug 2017
"""

import sys, os, subprocess, math

from pyned import find_differences, uncode, format_line
import read_config

#------------------------------------------------------------------------

import matplotlib.pyplot as plt


def find_indices(line):
    a = 1
    while line[a] != 'nan':
        a += 1
    n = a+1
    while line[n] != 'nan':
        n += 1
    return range(a+3,n), range(n+1, len(line))


def extract(data):
    R = data[1]   # radius
    M = data[4]   # microtubule polymer
    F = data[6]   # force
    B = data[8]   # number of times broken
    T = data[9]   # thickness
    return (R, M, B, T)


def make_plot(arg):
    fts = 14
    fig = plt.figure(figsize=(8, 6))
    # plot each data point with a symbol:
    for data in arg:
        (R, M, B, T) = extract(data)
        L = 2*math.pi*R
        S = 'x'
        if B < 3:
            S = '^'
        if B < 2:
            S = 's'
        if B == 0:
            plt.plot(L, M, 'o', markersize=L/2, markeredgewidth=T/10, markerfacecolor='white', markeredgecolor='black')
        else:
            plt.plot(L, M, S, markersize=4, markerfacecolor='none', markeredgecolor='blue')
    if 1:
        # add a R^4 curve going through (R=1, M=100)
        ref = ( 1, 100 )
        val = range(100)
        R = [ 0.5+0.02*x for x in val ]
        M = [ ref[1]*(x/ref[0])**4 for x in R ]
        L = [ 2*math.pi*x for x in R ]
        plt.plot(L, M, '-', color='blue', linewidth=1)    
    # label axes
    plt.ylabel('Available polymer length', fontsize=fts)
    plt.xlabel('Platelet perimeter (2 pi radius)', fontsize=fts)
    plt.xlim(3, 15)
    plt.ylim(45, 360)
    #plt.ylabel('force exterted by ring', fontsize=fts)
    plt.title('Ring simulations', fontsize=fts)
    fig.tight_layout()


#------------------------------------------------------------------------

def get_average(file, start):
    """
    Calculate mean of values beyond line 'start'
    """
    num = 0
    res = 0.0
    cnt = 0.0
    for line in file:
        num += 1
        if num > start:
            cnt += 1
            res += float(line)
    if cnt > 0:
        return res / cnt
    else:
        return 0


def get_ring(file, start):
    """
    Extract number of zero values beyond line 'start'
    """
    num = 0
    fail = 0
    mean = 0.0
    cnt = 0.0
    for line in file:
        num += 1
        if num > start:
            cnt += 1
            s = line.split()
            fail += ( int(s[0]) == 0 )
            mean += float(s[0])
    if cnt > 0:
        return ( fail, mean / cnt )
    else:
        return ( 0, 0 )


def get_parameters(path):
    pile = read_config.parse(path)
    res = {}
    try:
        cmd = read_config.get_command(pile, ['set', 'fiber', 'microtubule']);
        res = cmd.value('total_polymer')
    except:
        res = find_differences('config.cym', path+'/config.cym')
    return res
    
    
#------------------------------------------------------------------------

def process(path):
    """ 
        This extracts parameters from the config file,
        and values from files 'ring.txt' and 'force.txt'
    """
    if path.startswith('run'):
        res = [ path[3:] ]
    else:
        res = [ path ]
    if not os.path.isfile(path+'/config.cym'):
        return []        
    try:
        #par = get_parameters(path+'/config.cym')
        par = find_differences('config.cym', path+'/config.cym')
    except IOError as e:
        par = 'par_failed'
    res.extend(par)
    os.chdir(path)
    if not os.path.isfile('objects.cmo'):
        return []
    if not os.path.isfile('properties.cmo'):
        return []
    res.append('nan')
    # get forces:
    filename = 'force.txt'
    if not os.path.isfile(filename):
        subprocess.call(['report3', 'fiber:confinement', 'verbose=0'], stdout=open(filename, 'w'))
    with open(filename, 'r') as f:
        data = get_average(f, 0)
        res.append(data)
    res.append('nan')
    # get ring integrity data:
    filename = 'ring.txt'
    if not os.path.isfile(filename):
        subprocess.call(['report3', 'ring', 'verbose=0'], stdout=open(filename, 'w'))
    with open(filename, 'r') as f:
        data = get_ring(f, 0)
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
        sys.stderr.write("  Error: mising comparison base `config.cym'\n")
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
        plt.savefig('rings.pdf', dpi=300)
        #plt.show()
        plt.close()


#------------------------------------------------------------------------

if __name__ == "__main__":
    if len(sys.argv) < 2 or sys.argv[1].endswith("help"):
        print(__doc__)
    else:
        main(sys.argv[1:])

