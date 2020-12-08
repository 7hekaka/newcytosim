#!/usr/bin/env python3
#
# Make a plot using matplotlib
#
# F. Nedelec, 21 Aug 2017
#


"""
    Analysis script for the network contraction
    
Syntax:
    
    plot_bridges.py DIRECTORY_PATH
    
Description:
    
"""

import sys, os, math, subprocess, random, copy
import numpy as np
import matplotlib.pyplot as plt

from pyned import uncode
import read_config

#import matplotlib
#matplotlib.use('SVG')


fts = 14

#------------------------------------------------------------------------

def square(x):
    return x * x;

def get_meshsize(filename):
    """
        This extracts parameters from the config file, and compute the meshsize
    """
    pile = read_config.parse(filename)
    com = read_config.get_command(pile, ['set', 'space', '*'])
    geo = com.value("geometry")
    space_radius = float(geo.split()[1])
    space_volume = math.pi * square(space_radius)
    #print('space_radius  = ', space_radius)
    com = read_config.get_command(pile, ['new', 'fiber', '*'])
    fiber_length = com.value("length")
    nb_fiber = com.cnt
    p0 = 1.0 / square(math.pi) - 0.0235 * fiber_length / space_radius; # we previously used 0.09;
    nb_crossings = p0  * nb_fiber * ( nb_fiber - 1 ) * square( fiber_length / space_radius )
    mesh_size = fiber_length / ( 1 + 2 * nb_crossings / nb_fiber )
    #print('nb_fibers     = ', nb_fiber)
    #print('fibers_length = ', fiber_length)
    #print('nb_crossings  = ', nb_crossings)
    #print('mesh_size     = ', mesh_size)
    return mesh_size

#------------------------------------------------------------------------


def get_bridges(file):
    """
        Gather length and speed of bridges
    """
    dis = [];
    vel = [];
    for line in file:
        if len(line) < 4 or line[0] == '%':
            continue
        s = uncode(line).split()
        if len(s) == 4:
            dis.append(float(s[0]))
            vel.append(float(s[1]))
    return dis, vel


def plot_bridges(dis, vel):
    fig = plt.figure(figsize=(4, 2.5))
    plt.hist(dis, 32, facecolor='b')
    try:
        m = get_meshsize('config.cym')
        #print('meshsize=%.3f' % m)
        plt.plot([m, m], plt.ylim(), 'k-', linewidth=1)
    except Exception as e:
        print(e)
        pass
    plt.xlabel('Bridge length (um)', fontsize=fts)
    plt.ylabel('Histogram', fontsize=fts)
    #plt.title('Bridge length', fontsize=fts)
    fig.tight_layout()
    plt.savefig('bridges', dpi=150)
    #plt.show()
    plt.close()


def parse(dirpath):
    """
    Work in current directory
    """
    filename = 'bridges.txt'
    if not os.path.isfile(filename):
        subprocess.call(['reportN', 'network:bridges', 'frame=0'], stdout=open(filename, 'w'))
    with open(filename, 'r') as f:
        dis, vel = get_bridges(f)
    plot_bridges(dis, vel)


#------------------------------------------------------------------------

def main(args):
    paths = []
    for arg in args:
        if os.path.isdir(arg):
            paths.append(arg)
        else:
            sys.stderr.write("  Error: unexpected argument `%s'\n" % arg)
            sys.exit()
    if not paths:
        parse('.')
    else:
        cdir = os.getcwd()
        for p in paths:
            os.chdir(p)
            #sys.stdout.write('- '*32+p+"\n")
            parse(p)
            os.chdir(cdir)


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1].endswith("help"):
        print(__doc__)
    else:
        main(sys.argv[1:])

