#!/usr/bin/env python3
#
# Make a plot using matplotlib
#
# F. Nedelec, 19--21 February 2016
#


"""
    Analysis script for the network contraction
    
Syntax:
    
    plot_network.py DIRECTORY_PATH
    
Description:
    
"""

import sys, os, math, subprocess, random, copy
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import read_config

#matplotlib.use('SVG')


fts = 13
scale = 0.5

#------------------------------------------------------------------------

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


def differentiate(size):
    """
        Calculate derivative, 
        Returned array is one slot smaller
    """
    diff = []
    for i in range(1,len(size)):
        if size[i-1] > 0:
            diff.append((size[i]-size[i-1])/size[i-1])
        else:
            diff.append(0)
    return diff

#------------------------------------------------------------------------

def get_meshsize(filename):
    """
        This extracts parameters from the config file, and compute meshsize
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


def plot_scatter(size, cord):
    """
        Make a scatter-plot with data provided as arguments
    """
    xdata = np.array(size)
    ydata = np.array(cord)
    fig = plt.figure(figsize=(4, 4))
    ax = fig.add_axes([0.25, 0.17, 0.7, 0.7])
    plt.scatter(xdata, ydata)
    plt.plot([-scale, scale], [-scale, scale], 'g', linestyle='dotted')
    plt.plot([-scale, scale], [scale, -scale], 'g', linestyle='dotted')
    plt.xlim(-scale, scale)
    plt.ylim(-scale, scale)
    ax.set_xticks([-scale, 0, scale])
    ax.set_yticks([-scale, 0, scale])
    plt.title('Scatter plot', fontsize=fts)
    #fig.tight_layout()


def plot_intersections(tim, data):
    """
        Plot evolution in time of the bridge's length
    """
    fig = plt.figure(figsize=(5, 4))
    for k in data.keys():
        RGB = random_color()
        val = np.array(data[k])
        plt.plot(tim, val, '-', linewidth=0.5, color=RGB)
    plt.ylim(0, 5)
    plt.xlabel('Time (?)', fontsize=fts)
    plt.ylabel('Abscissa (um)', fontsize=fts)
    plt.title('Intersections', fontsize=fts)


def get_intersections(filename):
    """
        Calculate the absissa of bridges over time
    """
    tim = []
    res = []
    T  = 0
    with open(filename, 'r') as f:
        for line in f:
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
                for i in range(2, len(spl), 4):
                    a = float(spl[i])
                    g = int(spl[i+1])
                    k = szudzik(f,g)
                    if k in res:
                        res[k].append(a)
                    else:
                        res[k] = [a]
    return tim, res


def plot_bridge_count(tim, data):
    # get all keys:
    keys = set()
    for d in data:
        keys = keys.union(d.keys())
    #print(keys)
    fig = plt.figure(figsize=(5, 3))
    for k in keys:
        val = []
        for d in data:
            if k in d:
                val.append(d[k])
            else:
                val.append(0)
        #print(k, val)
        plt.plot(tim, val, '-', linewidth=1.5, color=random_color())
    #plt.ylim(0, max(yy))
    plt.xlabel('Time (s)', fontsize=fts)
    plt.ylabel('Nb of configurations', fontsize=fts)
    plt.title('Bridge counts', fontsize=fts)
    fig.tight_layout()
    plt.savefig('bridge_count', dpi=150)
    #plt.show()
    plt.close()


def plot_bridge_length(tim, val):
    fig = plt.figure(figsize=(5, 3))
    plt.plot(tim, val, 'b-', linewidth=3)
    try:
        m = get_meshsize('config.cym')
        plt.plot([min(tim), max(tim)], [m, m], 'k-', linewidth=1)
    except:
        pass
    plt.ylim(0, max(val)+min(val))
    plt.xlabel('Time (s)', fontsize=fts)
    plt.ylabel('distance between connectors', fontsize=fts)
    plt.title('Bridge length', fontsize=fts)
    fig.tight_layout()
    plt.savefig('bridge_length', dpi=150)
    #plt.show()
    plt.close()


def categorize(c):
    return c


def count_bridges(filename):
    """
        Count all the bridges
    """
    T = 0
    tim = []
    dis = []
    res = []
    cnt = {}
    T = 0
    with open(filename, 'r') as f:
        for line in f:
            spl = uncode(line).split()
            if len(spl) < 2:
                pass
            elif spl[0] == '%':
                if spl[1] == "start" or spl[1] == "time":
                    cnt = {}
                    s = 0
                    d = 0
                    T = float(spl[2])
                elif spl[1] == "end":
                    if s:
                        #print('number of bridges', s)
                        dis.append(d/s)
                        tim.append(T)
                        res.append(copy.deepcopy(cnt))
            elif len(spl) > 4:
                a = -1;
                for i in range(2, len(spl), 4):
                    b = a;
                    a = float(spl[i])
                    if b > 0:
                        d += a - b
                        s += 1
                    c = categorize(spl[i+2])
                    if c in cnt:
                        cnt[c] += 1
                    else:
                        cnt[c] = 1
    return tim, res, dis


def get_bridges(filename):
    """
        get the number of certain configurations found in bridges
    """
    res = {}
    res[0] = 0
    res[1] = 0
    res[2] = 0
    res[3] = 0
    res[4] = 0
    res[5] = 0
    with open(filename, 'r') as f:
        for line in f:
            spl = uncode(line).split()
            if len(spl) < 2:
                pass
            elif spl[0] == '%':
                if spl[1] == "end":
                    res[0] += 1
                    #print(res[1], res[2], res[3], res[4], res[5])
            elif len(spl) > 4:
                f = int(spl[1])
                for i in range(2, len(spl), 4):
                    a = float(spl[i])
                    g = int(spl[i+1])
                    c = spl[i+2]
                    if c in res:
                        res[c] += 1
                    else:
                        res[c] = 1
                    k = 5
                    if c == "0+1":
                        k = 1;
                    elif c.startswith("0+"):
                        k = 2;
                    elif c == "1+0":
                        k = 3
                    elif c.endswith("+0"):
                        k = 4
                    res[k] += 1
    return res


def triangle_surface(side):
    """
        Compute surface of a triangle using Hero's formula
    """
    P = 0.5 * sum(side)
    S2 = P*(P-side[0])*(P-side[1])*(P-side[2])
    if S2 > 0:
        return math.sqrt(S2)
    else:
        return 0


def get_triangle(filename):
    """
        Get parameters of a triangle from cytosim's report file
    """
    data = []
    side = [0, 0, 0]
    with open(filename, 'r') as f:
        for line in f:
            spl = uncode(line).split()
            if len(spl) < 2:
                pass
            elif spl[0] == '%':
                if spl[1] == "end":
                    data.append(triangle_surface(side))
            elif len(spl) > 7:
                f = int(spl[1])
                a = float(spl[2])
                b = float(spl[6])
                if 0 < f and f < 4:
                    side[f-1] = b - a
    return data


def get_cords(filename):
    """
        extract lengths of bridges
    """
    T = 0
    time = []
    acum = []
    dcum = []
    past = {}
    cabs = {}
    data = []
    ncor = []
    pcor = []
    with open(filename, 'r') as f:
        for line in f:
            spl = uncode(line).split()
            if len(spl) < 2:
                pass
            elif spl[0] == '%':
                if spl[1] == "start" or spl[1] == "time":
                    T = float(spl[2])
                    past = cabs
                    cabs = {}
                    acum = []
                    dcum = []
                elif spl[1] == "end":
                    n = len(acum)
                    if n:
                        time.append(T)
                        #print(dcum)
                        ncum = [i for i in dcum if i < 0]
                        a = sum(acum)/n
                        da = sum(dcum)/n
                        na = sum(ncum)/n
                        data.append(da/a)
                        ncor.append(na/a)
            elif len(spl) > 6:
                f = int(spl[1])
                for i in range(6, len(spl), 4):
                    a = float(spl[i-4])
                    g = int(spl[i-3])
                    b = float(spl[i])
                    h = int(spl[i+1])
                    cabs[szudzik(f, g)] = a;
                    cabs[szudzik(f, h)] = b;
                    da = b - a
                    #print(f, g, h, a, b)
                    try:
                        oa = past[szudzik(f, h)] - past[szudzik(f, g)]
                        acum.append(da)
                        dcum.append(da-oa)
                    #print(f, g, h, round(oa,6), round(da,6))
                    except:
                        pass
    return [ time, data, ncor ]


def plot_size(time, data):
    """
        Plot size as a function of time
    """
    fig = plt.figure(figsize=(5, 4))
    plt.plot(time, data, 'b-', linewidth=3)
    s0 = data[0]
    plt.plot([min(time), max(time)], [s0, s0], 'k-', linewidth=0.5)
    #plt.xlim(0, 100)
    plt.ylim(0, min(data)+max(data))
    plt.xlabel('Time (s)', fontsize=fts)
    plt.ylabel('Radius (um)', fontsize=fts)
    plt.title('Network size', fontsize=fts)
    fig.tight_layout()


def get_moment(file):
    """
        Get network size as a function of time
        """
    tim = []
    mom = []
    T = 0
    M = 0
    for line in file:
        s = uncode(line).split()
        if len(s) < 2:
            pass
        elif s[0] == '%':
            if s[1] == "start" or s[1] == "time":
                T = float(s[2])
            elif spl[1] == "moment":
                M = float(s[-1])
            elif s[1] == "end":
                tim.append(T)
                mom.append(M)
        elif len(s) == 9 and s[0].isalpha():
            M = float(s[8])
        elif len(s) == 10 and s[1].isalpha():
            T = []
            M = []
            tim.append(float(s[0]))
            mom.append(float(s[9]))
    return tim, mom


def get_size(file):
    tim, mom = get_moment(file)
    siz = [ math.sqrt(2*m) for m in mom ]
    return tim, siz


def parse(dirpath):
    """
    Work in current directory
    """
    filename = 'connectors.txt'
    if not os.path.isfile(filename):
        subprocess.call(['reportN', 'fiber:connector'], stdout=open(filename, 'w'))
    if 1:
        time, data, dis = count_bridges(filename)
        plot_bridge_count(time, data)
        plot_bridge_length(time, dis)
    if 0:
        time, cord, ncord = get_cords(filename)
        with open(filename, 'r') as f:
            time, data = get_size(f)
    if 0:
        plot_scatter(differentiate(size), cord)
        plt.plot(np.array(differentiate(size)), np.array(ncord), 'r.')
        plt.title('Radius from connectors', fontsize=fts)
        plt.xlabel('dR/R', fontsize=fts)
        plt.ylabel('da/a', fontsize=fts)
        plt.savefig('scatter_connect', dpi=150)
        plt.close()
    if 1:
        with open(filename, 'r') as f:
            time, data = get_size(f)
        prune_values(time, size)
        plot_size(time, size)
        plt.title('Radius from connectors', fontsize=fts)
        plt.savefig('size', dpi=150)
        plt.close()
    if 0:
        tim, data = get_intersections(filename)
        plot_intersections(tim, data)
        plt.savefig('connectors_abscissa', dpi=150)
#        plt.show()
        plt.close()
    if 0:
        data = get_connectors(filename)
        print(data[1], data[2], data[3], data[4], data[5])
    if 0:
        surf = get_triangle(filename)
        plot_scatter(differentiate(surf), cord)
        plt.xlabel('d(triangle_surface)/triangle_surface', fontsize=fts)
        plt.ylabel('da/a', fontsize=fts)
        plt.savefig('triangle', dpi=150)
#        plt.show()
        plt.close()
    if 0:
        filename = 'cross.txt'
        if not os.path.isfile(filename):
            subprocess.call(['report2', 'fiber:intersection', 'details=0'], stdout=open(filename, 'w'))
        with open(filename, 'r') as f:
            time, data = get_size(f)
        print(differentiate(size))
        plot_scatter(differentiate(size), cord)
        plt.title('Radius from crossings', fontsize=fts)
        plt.xlabel('dR/R', fontsize=fts)
        plt.ylabel('da/a', fontsize=fts)
        plt.savefig('scatter_cross', dpi=150)
        plt.close()
    if 0:
        plot_size(time, size)
        plt.title('Radius from crossings', fontsize=fts)
        plt.savefig('size_cross', dpi=150)
        plt.close()
    if 0:
        filename = 'mom.txt'
        if not os.path.isfile(filename):
            subprocess.call(['report2', 'fiber:moment'], stdout=open(filename, 'w'))
        with open(filename, 'r') as f:
            time, data = get_size(f)
        print(differentiate(size))
        plot_scatter(differentiate(size), cord)
        plt.title('Radius from points', fontsize=fts)
        plt.xlabel('dR/R', fontsize=fts)
        plt.ylabel('da/a', fontsize=fts)
        plt.savefig('scatter_mass', dpi=150)
#       plt.show()
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

