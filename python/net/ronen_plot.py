#!/usr/bin/env python3
#
# A script to plot for project with Ronen Zaidel-Bar
#
# F. Nedelec, Strasbourg, 16.12.2021, 8-12.1.2022, 14--16.3.2022


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

def modifs(mod):
    keys = ['Reference', ' 2xActin', ' 1.33xArp23', ' 2xMyosin', 'Dead Arp23', '+Capping', '++Capping', '33% Reduced Actin']
    res = ''
    if not mod:
        return keys[0]
    if mod == 99:
        return '2xActin + 2xMyosin'
    for i in { 0, 1, 2, 3, 4, 5, 6 }:
        if int((mod>>i)&1):
            res+=keys[i+1]
    return res


def one_plot(data, X, Y, dd = 0):
    pX = [ d[0].strip("run") for d in data if d[4]==X ]
    pY = [ d[0].strip("run") for d in data if d[4]==Y ]
    dX = [ d[2] for d in data if d[4]==X ]
    dY = [ d[2] for d in data if d[4]==Y ]
    if type(dd) is list:
        dY = dd
    fig = plt.figure(figsize=(3.84, 3.84)) #(7,7))
    plt.scatter(dX, dY, marker='o', s=8, c='blue')
    # add diagonal:
    M = math.ceil(20*max(dX+dY))*0.05
    plt.plot([0, M], [0, M], 'k:', linewidth=1)
    plt.xlim(0, M)
    plt.ylim(0, M)
    if type(dd) is int and dd > 0:
        for i in range(len(dX)):
            t = pX[i]+','+pY[i]
            plt.text(dX[i], dY[i], t, size=5)
    plt.xlabel('%s' % modifs(X), fontsize=fts)
    plt.ylabel('%s' % modifs(Y), fontsize=fts)
    plt.title('Contraction rate (um/s)', fontsize=fts)
    fig.tight_layout()
    plt.savefig('0_contraction%i%i.png' %(X,Y), dpi=500)
    plt.close()


def two_plot(data, X, Y):
    pL = [ d[5] for d in data if d[4]==X ]
    dX = [ d[2] for d in data if d[4]==X ]
    dY = [ d[2] for d in data if d[4]==Y ]
    fig = plt.figure(figsize=(3.84, 3.84)) #(7,7))
    for i in range(len(dX)):
        plt.scatter(dX[i], dY[i], marker='o', s=pL[i]/256, c='#0000FF', edgecolors='#AAAAFF', linewidths=0.25)
    # add diagonal:
    M = math.ceil(20*max(dX+dY))*0.05
    plt.plot([0, M], [0, M], 'k:', linewidth=1)
    plt.xlim(0, M)
    plt.ylim(0, M)
    plt.xlabel('%s' % modifs(X), fontsize=fts)
    plt.ylabel('%s' % modifs(Y), fontsize=fts)
    plt.title('Contraction rate (um/s)', fontsize=fts)
    fig.tight_layout()
    plt.savefig('2_contraction%i%i.png' %(X,Y), dpi=500)
    plt.close()


def plot_flen(data, X, Y):
    dX = [ d[5] for d in data if d[4]==X ]
    dY = [ d[5] for d in data if d[4]==Y ]
    fig = plt.figure(figsize=(3.84, 3.84))
    plt.scatter(dX, dY, marker='o', s=8, c='blue')
    # add diagonal:
    M = math.ceil(20*max(dX+dY))*0.05
    plt.plot([0, M], [0, M], 'k:', linewidth=1)
    plt.xlim(0, M)
    plt.ylim(0, M)
    plt.xlabel('%s' % modifs(X), fontsize=fts)
    plt.ylabel('%s' % modifs(Y), fontsize=fts)
    plt.title('Polymer mass (um)', fontsize=fts)
    fig.tight_layout()
    plt.savefig('1_polymer%i%i.png' %(X,Y), dpi=500)
    plt.close()


def many_plots(data):
    """
        Make summary plots
    """
    #P, A, B, C, M, L = zip(*data)
    # get unique values:
    pool = [ d[4] for d in data ]
    unique = set(pool)
    for m in unique:
        print(m, ': ', pool.count(m))
    # check all pairwise combinations:
    for X in { 0 }:
        for Y in unique:
            if Y > X:
                one_plot(data, X, Y)
                plot_flen(data, X, Y)
    two_plot(data, 0, 64)
    # check some pairs:
    if { 5, 7 }.issubset(unique):
        one_plot(data, 5, 7)
    # check additivity of the effect:
    if { 1, 4, 5 }.issubset(unique):
        # check combined effects of 1 and 4 and compare against 5:
        c0 = [ d[2] for d in data if d[4]==0 ]
        c1 = [ d[2] for d in data if d[4]==1 ]
        c4 = [ d[2] for d in data if d[4]==4 ]
        pp = [ b * c / a for a, b, c in zip(c0, c1, c4) ]
        one_plot(data, 5, 99, pp)

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
            data.append([s[0], float(s[1]), float(s[2]), float(s[3]), int(s[4]), float(s[5])])
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
        files.append(p+'/rates.txt')
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

