#!/usr/bin/env python3
#
# A script to plot for project with Ronen Zaidel-Bar
#
# F. Nedelec, Strasbourg, 16.12.2021, 8-12.1.2022


"""
    Read all contraction rates in 'rates.txt' and make master plot
    
Syntax:
    
    ronen_meta.py DATA_FILE
    
Description:
    
"""

plot_ratio = 0
condition = 0
#font size:
fts = 10

import sys, os, math
try:
    import matplotlib
    #matplotlib.use('SVG')
    import matplotlib.pyplot as plt
except:
    print("  Error: could not load matplotlib")
    sys.exit()


def get_parameters(filename):
    """
        Read from config file numeric data associated with 'new'
    """
    res = {}
    with open(filename, 'r') as f:
        for line in f:
            s = line.split()
            if len(s) == 3:
                if s[0] == 'total_polymer':
                    res[s[0]] = int(s[2])
                elif s[0] == 'new' and s[1].isdigit():
                    res[s[2]] = int(s[1])
    return res

#-------------------------------------------------------------------------------

def normalize(L):
    print("%i  Z values in [ %f %f ]" % (len(L), min(L), max(L)), end='')
    print(" %i negative" % len([ 1 for x in L if x < 0 ]), end='')
    print(" %i positive" % len([ 1 for x in L if x > 0 ]), end='')
    #print(sorted([round(x,2) for x in L]))
    a = 1.0 / max(L)
    T = [ x * a for x in L ]
    print(" normalized to [ %f %f ]" % (min(T), max(T)))
    return T

def transform_value(x, a, b):
    if x > 0:
        return x * a
    else:
        return x * b

def transform(Z):
    L = [ math.log(max(x, 0.01)) for x in Z ]
    print("%i  log(Z) values in [ %f %f ]" % (len(L), min(L), max(L)), end='')
    print(" %i negative" % len([ 1 for x in L if x < 0 ]), end='')
    print(" %i positive" % len([ 1 for x in L if x > 0 ]), end='')
    #print(sorted([round(x,2) for x in L]))
    a = 1.0 # / max(L)
    b = 2.0 # / min(L)
    T = [ transform_value(x, a, b) for x in L ]
    print(" transformed to [ %f %f ]" % (min(T), max(T)))
    return T

def cap(data):
    res = [ min(x, 1) for x in data ]
    res = [ max(x, -1) for x in res ]
    return res


def split(X, Y, Z, i):
    "split data according to sign of i-th value"
    xyz = list(zip(X, Y, Z))
    N = [ x for x in xyz if x[i] < 0 ]
    if N:
        N = list(zip(*N))
        if i > 0:
            N[i] = [ -x for x in N[i] ]
        P = [ x for x in xyz if x[i] >= 0 ]
        P = list(zip(*P))
    else:
        P = (X, Y, Z)
    return N, P


#-------------------------------------------------------------------------------

def plot(X, Y, Z):
    """
        plot data at position ( X, Y ) with color according to Z
    """
    fig = plt.figure(figsize=(5, 4))
    N, P = split(X, Y, Z, 2)
    if N:
        plt.scatter(N[0], N[1], marker='s', s=12, c=N[2], cmap='Blues', label='Less contraction')
    if P:
        plt.scatter(P[0], P[1], marker='o', s=16, c=P[2], cmap='Oranges', label='More contraction')
    # get limits
    mX = math.ceil(max(X)/10) * 10
    mY = math.ceil(max(Y)/10) * 10
    #M = max(mX, mY)
    #plt.plot([0, M], [0, M], 'k-', linewidth=1)
    plt.xlim(0, mX)
    plt.ylim(0, mY)
    return fig


def one_plot(params, i, j, values):
    #print(params, i, j)
    X = params[i]
    Y = params[j]
    fig = plot(X, Y, values)
    plt.xlabel(i, fontsize=fts)
    plt.ylabel(j, fontsize=fts)
    if plot_ratio:
        plt.title('Condition %i / condition 0' % condition, fontsize=fts)
    else:
        plt.title('Condition %i' % condition, fontsize=fts)
    plt.legend(loc='upper right', fontsize=7)
    fig.tight_layout()
    if plot_ratio:
        plt.savefig('2_%s_%s.png' %(i,j), dpi=150)
    else:
        plt.savefig('1_%s_%s.png' %(i,j), dpi=150)
    plt.close()


def many_plots(data):
    """
        Make summary plots
    """
    global condition
    P, A, B, C, M = zip(*data)
    # get unique values:
    pool = {}
    for m in set(M):
        pool[m] = [ d for d in data if d[4]==m ]
        print(m, ': ', len(pool[m]))
    # fix error in one screen by reducing to smallest common size:
    if len(pool[2]) == 2 * len(pool[0]):
        pool[2] = pool[2][0::2]
    keys = get_parameters(P[0]+'/config.cym').keys()
    print(list(keys))
    params = {}
    for k in keys:
        params[k] = []
    values = []
    # adjust index of condition to consider:
    if not condition in pool:
        condition = 1
    # load parameters for all datasets:
    if plot_ratio:
        for p0, p1 in zip(pool[0], pool[condition]):
            values.append(p1[2]/p0[2]) # ratio of the rates
            pam = get_parameters(p0[0]+'/config.cym')
            #pim = get_parameters(p1[0]+'/config.cym')
            for k, v in pam.items():
                params[k].append(v)
        values = cap(transform(values))
    else:
        for p0 in pool[condition]:
            values.append(p0[2]) # rate of reference system
            pam = get_parameters(p0[0]+'/config.cym')
            
            for k, v in pam.items():
                params[k].append(v)
        values = normalize(values)
    # check pairs:
    for i in keys:
        for j in keys:
            if j > i:
                one_plot(params,i,j,values)

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

