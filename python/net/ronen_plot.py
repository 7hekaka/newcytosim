#!/usr/bin/env python3
#
# A script to plot for project with Ronen Zaidel-Bar
#
# F. Nedelec, Strasbourg, 16.12.2021, 8-9.1.2022


"""
    Read all contraction rates in 'scores.txt' and make master plot
    
Syntax:
    
    ronen_plot.py DATA_FILE
    
Description:
    
"""

order = 7
#font size:
fts = 10

import sys, os, math
try:
    import matplotlib
    #matplotlib.use('SVG')
    import matplotlib.pyplot as plt
except:
    print("  Error: could not load matplotlib in python " + sys.version)
    sys.exit()

#-------------------------------------------------------------------------------

def read_data(filename):
    """
        Read numeric data from file
    """
    data = []
    with open(filename, 'r') as f:
        for line in f:
            if len(line) > 5:
                s = line.split()
                if s[0] != '%':
                    data.append([s[0], float(s[1]), float(s[2]), float(s[3]), int(s[4])])
    return data


def plot(poolX, poolY, modX, modY, title='Rate correlation'):
    """
        Make one plot to compare data in conditions X and Y
    """
    fig = plt.figure(figsize=(5, 4))
    plt.plot(poolX, poolY, 'bo', markersize=3)
    # add diagonal:
    M = math.ceil(20*max(poolX))*0.05
    plt.plot([0, M], [0, M], 'k-', linewidth=1)
    plt.xlim(0, M)
    plt.ylim(0, M)
    plt.xlabel('Contraction (%s)' % modX, fontsize=fts)
    plt.ylabel('Contraction (%s)' % modY, fontsize=fts)
    plt.title(title, fontsize=fts)
    fig.tight_layout()


def modifs(mod):
    keys = ['Reference', 'MoreActin', 'LessArp23', 'MoreMyosin']
    res = ''
    if not mod:
        return keys[0]
    mod1=int(mod&1)
    mod2=int((mod>>1)&1)
    mod4=int((mod>>2)&1)
    if mod1:
        res+=keys[1]
    if mod2:
        res+=keys[2]
    if mod4:
        res+=keys[3]
    return res


def one_plot(pool, X, Y):
    plot(pool[X], pool[Y], modifs(X), modifs(Y));
    plt.savefig('result%i%i.png' %(X,Y), dpi=150)
    plt.close()


def many_plots(data, ord):
    """
        Make summary plots
    """
    print("Master plot with %i datapoints" % len(data))
    P, A, B, C, M = zip(*data)
    # get unique values:
    pool = {}
    for m in set(M):
        P = [ d[2] for d in data if d[4]==m ]
        print(m, ': ', len(P))
        pool[m] = P
    
    for X in { 0 }:
        for Y in { 1, 2, 4, 5, 7 }:
            one_plot(pool, X, Y)
    one_plot(pool, 5, 7)

#-------------------------------------------------------------------------------

def main(args):
    paths = []
    files = []
    for arg in args:
        if os.path.isdir(arg):
            paths.append(arg)
        elif os.path.isfile(arg):
            files.append(arg)
        else:
            sys.stderr.write("  Error: unexpected argument `%s'\n" % arg)
            sys.exit()
    for p in paths:
        files.append(p+'/scores.txt');
    if not files:
        files = ['scores.txt']
    for f in files:
        data = read_data(f)
        many_plots(data, order)


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1].endswith("help"):
        print(__doc__)
    else:
        main(sys.argv[1:])

