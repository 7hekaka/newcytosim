#!/usr/bin/env python3
#
# get_network.py
#
# simple extraction of data from files
#
# Copyright F. Nedelec, 24.08.2017

"""
    Visit given run directories, and predict contraction from bridges characteristics

Syntax:

    get_network.py DIRECTORIES
    
Example:

    get_network run???? > net.txt

Description:

    This script analyzes simulations of network contraction.


F. Nedelec, 21.08.2017
"""

import sys, os, subprocess, math
from pyned import find_differences, uncode, format_line

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


def make_plot(arg):
    fts = 14
    fig = plt.figure(figsize=(7, 5))
    xinx = 1
    xval = []
    yval = []
    inxT, inxD = find_indices(arg[0])
    arg = sorted(arg, key=lambda line:line[xinx])
    # transpose list of lists:
    data = list(zip(*arg))
    # plot the prediction with line:
    for i in inxT:
        c = inxT.index(i) * 0.75 / len(inxT)
        plt.plot(data[xinx], data[i], 'k-', linewidth=2, color=(c,c,c))
    # plot all the contraction data points:
    for i in inxD:
        c = inxD.index(i) * 0.75 / len(inxD)
        plt.plot(data[xinx], data[i], 'o', markersize=5, markerfacecolor=(1,1,1,0), markeredgecolor=(c,c,1-c))
    # label axes
    plt.xlabel('Column 1', fontsize=fts)
    plt.ylabel('Relative Contraction Rate (s-1)', fontsize=fts)
    plt.title('Bridges vs. Contraction', fontsize=fts)
    fig.tight_layout()
    plt.savefig('net', dpi=150)
    #plt.show()
    plt.close()


#------------------------------------------------------------------------

def count_bridges(file):
    """
        Average length and speed of bridges
    """
    cnt = 0;
    dis = 0;
    vel = 0;
    for line in file:
        if len(line) < 4 or line[0] == '%':
            continue
        s = uncode(line).split()
        if len(s) == 4:
            d = float(s[0])
            v = float(s[1])
            if v > 0:
                v = 0;
            cnt += 1
            dis += d
            vel += v
    if cnt > 0:
        return ( dis / cnt, vel / cnt, vel / dis )
    else:
        return ( 0, 0, 0 )


def count_bridges_binned(file):
    """
        Average length and speed of bridges
    """
    #bins = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 1000]
    bins = [1.0, 2.0, 3.0, 4.0, 1000]
    nbin = len(bins)
    dis = [0] * nbin
    vel = [0] * nbin
    cnt = [0] * nbin
    for line in file:
        if len(line) < 4 or line[0] == '%':
            continue
        s = uncode(line).split()
        if len(s) == 4:
            d = float(s[0])
            v = float(s[1])
            test = [ 1 if d < x else 0 for x in bins ]
            try:
                inx = test.index(1)
            except:
                inx = nbin-1
            if v < 0:
                vel[inx] += v
            cnt[inx] += 1
            dis[inx] += d
    res = [0] * (nbin+2)
    for i, v in enumerate(vel):
        if cnt[i] > 0:
            res[i+2] = v / dis[i]
    return res


def get_contraction(file):
    """
        Extract relative contraction rate from 'mom.txt' file
    """
    T0 = 0
    Z0 = 0
    res = []
    for line in file:
        line = uncode(line)
        if len(line) < 4 or line[0] == '%':
            continue
        s = line.split()
        T = float(s[0])
        #Z = math.pi*2*float(s[-1])
        Z = math.sqrt(2*float(s[-1]))
        if T0 > 0:
            dZ = ( Z - Z0 ) / ( T - T0 )
            res.append(dZ/Z0)
        else:
            T0 = T
            Z0 = Z
    return res


#------------------------------------------------------------------------

def process(path):
    """ 
        This extracts parameters from the config file,
        and values from 'mom.txt'
    """
    # add directory path
    if path.startswith('run'):
        res = [ path[3:] ]
    else:
        res = [ path ]
    # get parameter values
    if not os.path.isfile(path+'/config.cym'):
        return []
    try:
        par = find_differences('config.cym', path+'/config.cym')
    except IOError as e:
        par = 'par_failed'
    res.extend(par)
    res.append('nan')
    os.chdir(path)
    if not os.path.isfile('objects.cmo'):
        return []
    if not os.path.isfile('properties.cmo') and not os.path.isfile('properties.cmp'):
        return []

    # get bridge characteristics:
    filename = 'bridges.txt'
    if not os.path.isfile(filename):
        subprocess.call(['reportN', 'network:bridges', 'frame=0'], stdout=open(filename, 'w'))
    with open(filename, 'r') as f:
        data = count_bridges_binned(f)
        res.extend(data)
    res.append('nan')

    # get contraction rate from file 'mom.txt'
    filename = 'mom.txt'
    if not os.path.isfile(filename):
        subprocess.call(['reportN', 'fiber:moment', 'prefix=time'], stdout=open(filename, 'w'))
    with open(filename, 'r') as f:
        data = get_contraction(f)
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


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1].endswith("help"):
        print(__doc__)
    else:
        main(sys.argv[1:])


