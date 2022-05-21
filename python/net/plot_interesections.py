#!/usr/bin/env python3
#
# Make a plot using matplotlib
#
# F. Nedelec, 12 May 2016
#


"""
    Analysis script for the network contraction
    
Syntax:
    
    plot_intersections.py DIRECTORY_PATH
    
Description:
    
"""

import sys, os, math, random
import numpy as np
import matplotlib.pyplot as plt
import read_config
#matplotlib.use('SVG')


fts = 13
scale = 0.5

#------------------------------------------------------------------------

def square(x):
    return x * x;

def uncode(arg):
    try:
        if isinstance(arg, unicode):
            return str(arg.decode('utf-8'))
    except:
        pass
    return arg


def prune_values(time, data):
    """
        """
    i = len(data)-1
    while i >= 0:
        if not math.isfinite(data[i]):
            data.pop(i)
            time.pop(i)
        i-=1


def random_color():
    R=random.random()
    G=random.random()
    B=random.random()
    return (R, G, B)


def szudzik(a, b):
    """
        Matthew Szudzik's pairing function
    """
    if a >= b:
        return a * a + a + b
    else:
        return a + b * b


def kizduzs(z):
    """
        Reverse Matthew Szudzik's pairing function
    """
    b = math.floor(math.sqrt(z))
    a = z - b * b;
    if a < b:
        return [ a, b ]
    else:
        return [ b, a - b ]


def differentiate(data):
    """
        Calculate differences of successive data, 
        Returned array is one slot smaller
    """
    res = []
    for i in range(1,len(data)):
        res.append(data[i]-data[i-1])
    return res


def histogram(data, bin_size):
    """
        calculate histogram with regular bins of size bin_size
    """
    zer = min(data)
    val = [ math.floor(( x - zer ) / bin_size ) for x in data ]
    #print(max(val))
    res = [ val.count(n) for n in range(0, max(val)+1) ]
    pos = [ (n+0.5) * bin_size for n in range(0, max(val)+1) ]
    return pos, res

#------------------------------------------------------------------------

def get_meshsize(filename):
    """
        This extracts parameters from the config file, and computes the meshsize
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


def plot_histogram(data, counts):
    fig = plt.figure(figsize=(5, 4))
    plt.plot(data, counts, '.', linewidth=2, color=random_color())
    #plt.ylim(0, 5)
    plt.xlabel('Distance (um)', fontsize=fts)
    plt.ylabel('Counts', fontsize=fts)
    plt.title('Histogram', fontsize=fts)
    fig.tight_layout()


def get_intersections(filename):
    """
        Calculate an histogram of the distances between interesections
    """
    tim = []
    res = dict()
    T  = 0
    with open(filename, 'r') as f:
        for line in f:
            #print(line, end='')
            spl = uncode(line).split()
            if len(spl) < 2:
                pass
            elif spl[0] == '%':
                if spl[1] == "start" or spl[1] == "time":
                    T = float(spl[2])
                elif spl[1] == "end":
                    tim.append(T)
            elif len(spl) > 4:
                f = int(spl[1])
                a = float(spl[2])
                g = int(spl[3])
                b = float(spl[4])
                if f in res:
                    res[f].append(a)
                else:
                    res[f] = [a]
                if g in res:
                    res[g].append(b)
                else:
                    res[g] = [b]
    # sort all results:
    for f in res.keys():
        res[f].sort()
    return tim, res


def parse(dirpath):
    """
    Work in current directory
    """
    if 1:
        mesh_size = get_meshsize('config.cym')
        print('mesh_size     = ', mesh_size)
    filename = 'intersections.txt'
    if not os.path.isfile(filename):
        subprocess.call(['report', 'fiber:intersections', 'details=2'], stdout=open(filename, 'w'))
    tim, res = get_intersections(filename)
    #calculate the distance between consecutive interactions:
    data = []
    for f in res.keys():
         data.extend(differentiate(res[f]))
    print('mean_distance = ', sum(data)/len(data))
    #print(data)
    pos, counts = histogram(data, 0.002)
    plot_histogram(pos, counts)
    plt.plot([mesh_size, mesh_size], plt.ylim(), '-', linewidth=1)
    
    plt.title('Distance between intersections', fontsize=fts)
    plt.savefig('distances', dpi=150)
#   plt.show()
    plt.close()


#------------------------------------------------------------------------

def main(args):
    global scale
    paths = []
    
    for arg in args:
        if os.path.isdir(arg):
            paths.append(arg)
        elif arg.startswith('scale='):
            scale=float(arg[6:])
        else:
            sys.stderr.write("  Error: unexpected argument `%s'\n" % arg)
            sys.exit()
    
    if not paths:
        parse('.')
    else:
        cdir = os.getcwd()
        for p in paths:
            os.chdir(p)
            sys.stdout.write('- '*32+p+"\n")
            try:
                parse(p)
            except Exception as e:
                sys.stderr.write("Error: %s\n" % repr(e))
            os.chdir(cdir)


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1].endswith("help"):
        print(__doc__)
    else:
        main(sys.argv[1:])

