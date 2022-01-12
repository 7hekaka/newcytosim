#!/usr/bin/env python3
#
# A script to plot for project with Ronen Zaidel-Bar
#
# F. Nedelec, Strasbourg, 16.12.2021, 8-12.1.2022


"""
    Read all contraction rates in 'rates.txt' and make master plot
    
Syntax:
    
    ronen_plot.py DATA_FILE
    
Description:
    
"""

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

def plot(X, Y):
    """
        Make one plot to compare data in conditions X and Y
    """
    fig = plt.figure(figsize=(5, 4))
    plt.scatter(X, Y, marker='o', s=8, c='blue')
    # add diagonal:
    M = math.ceil(20*max(X))*0.05
    plt.plot([0, M], [0, M], 'k-', linewidth=1)
    plt.xlim(0, M)
    plt.ylim(0, M)
    return fig
    

def modifs(mod):
    keys = ['Reference', 'MoreActin', 'MoreArp23', 'MoreMyosin']
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
    fig = plot(pool[X], pool[Y]);
    plt.xlabel('Contraction (%s)' % modifs(X), fontsize=fts)
    plt.ylabel('Contraction (%s)' % modifs(Y), fontsize=fts)
    plt.title('Rate correlation', fontsize=fts)
    fig.tight_layout()
    plt.savefig('0_contraction%i%i.png' %(X,Y), dpi=150)
    plt.close()


def many_plots(data):
    """
        Make summary plots
    """
    P, A, B, C, M = zip(*data)
    # get unique values:
    pool = {}
    for m in set(M):
        pool[m] = [ d[2] for d in data if d[4]==m ]
        print(m, ': ', len(pool[m]))
    # reduce to smallest common size:
    if 2 in pool and len(pool[2]) == 2 * len(pool[0]):
        pool[2] = pool[2][0::2]
    # check all conditions:
    for X in { 0 }:
        for Y in pool.keys():
            if Y > X:
                one_plot(pool, X, Y)
    # check some pairs:
    if 7 in pool:
        one_plot(pool, 5, 7)
    if { 1, 4, 7 }.issubset(pool):
        #check combination of 1 and 4 against 5:
        pool[21] = [ b * c / a for a, b, c in zip(pool[0], pool[1], pool[4]) ]
        one_plot(pool, 5, 21)

#-------------------------------------------------------------------------------

def read_data(filename):
    """
        Read numeric data from file
    """
    data = []
    f = open(filename, 'r')
    for line in f:
        s = line.split()
        if len(line) > 4 and s[0] != '%':
            data.append([s[0], float(s[1]), float(s[2]), float(s[3]), int(s[4])])
    f.close()
    print("Collected %i datapoints" % len(data))
    return data


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
        files.append(p+'/rates.txt');
    if not files:
        files = ['rates.txt']
    for f in files:
        try:
            data = read_data(f)
            many_plots(data)
        except FileNotFoundError as e:
            sys.stderr.write(str(e)+'\n')


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1].endswith("help"):
        print(__doc__)
    else:
        main(sys.argv[1:])

