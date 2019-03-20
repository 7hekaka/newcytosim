#!/usr/bin/env python3
#
# get_platelet.py
#
# Copyright F. Nedelec, 2011 - 14.11.2018

"""
    Collect data from multiple given run directories, and print output to standard output
    Make several plots to summarize data as a function of platelet size

Syntax:

    get_platelet.py DIRECTORIES
    
Example:

    get_platelet.py run???? > scaling.txt

Description:

    This script is used to analyze platelet simulations

F. Nedelec, July 2017 -- 09.2018 -- 11.2018
"""

import sys, os, subprocess, math, copy, random
import matplotlib.pyplot as plt
from pyned import uncode, format_line, simple_linear_fit, linear_fit, power_fit, powerlaw_fit, exponential_fit, frange

title = os.path.basename(os.getcwd())

#------------------------------------------------------------------------
def cntminmax(arg):
    data = copy.deepcopy(arg)
    cnt = 0
    ix = math.inf
    sx = -math.inf
    for x, y in data:
        ix = min(ix, x)
        sx = max(sx, x)
        ++cnt
    return ( cnt, ix, sx )


def reservoir(arg, cnt):
    """ Reservoir sampling algorithm """
    i = 0
    res = []
    for i, x in enumerate(arg):
        if i < cnt:
            res[i] = x
        else:
            j = random.randint(0, i)
            if j <= cnt:
                res[j] = x
    return res


def bootstrap(arg):
    cnt, ix, sx = cntminmax(arg)
    data = list(copy.deepcopy(arg))
    for i in range(128):
        boot = [ random.choice(data) for x in data ]
        a, b = simple_linear_fit(boot)
        iv = a + b * ix
        sv = a + b * sx
        plt.plot((ix, sx), (iv, sv), '-', linewidth=0.5, color=[0,1,0], alpha=0.1)


def plot_data(data, power_fit, exp_fit):
    cnt, ix, sx = cntminmax(data)
    bootstrap(data)
    for x, y in data:
        plt.plot(x, y, 'o', markersize=4, markerfacecolor='none', markeredgewidth=0.5, markeredgecolor=[0,0,1], alpha=0.7)
    # add fitted power:
    bot, top = plt.ylim()
    if power_fit:
        scale, power = power_fit
        X = frange(ix, sx, 100)
        Y = [ scale*(x**power) for x in X ]
        plt.plot(X, Y, '-', color='blue', linewidth=0.5, alpha=0.5)
        plt.ylim(bot, top)
        plt.text(4, top, "fit %.3f R^%.3f\n" % (scale, power), fontsize=6, color='blue')
    if exp_fit:
        a, b = exp_fit
        X = frange(ix, sx, 100)
        Y = [ a*math.exp(b*x) for x in X ]
        plt.plot(X, Y, '-', color='red', linewidth=0.5, alpha=0.5)
        plt.ylim(bot, top)
        plt.text(15, top, "fit %.3f exp(%.3f x)\n" % (a, b), fontsize=6, color='red')


def laplace(data):
    fts = 14
    fig = plt.figure(figsize=(9, 6))
    for L, x, y in data:
        plt.plot(x, y, 'o', markersize=L/2, markerfacecolor='none', markeredgecolor='blue')
    plt.xlabel("ring tension", fontsize=fts)
    plt.ylabel("radial pressure", fontsize=fts)
    plt.title(title, fontsize=18)
    plt.savefig('p-laplace.pdf', dpi=300)
    plt.close()


def plot(data, label, power):
    fig = plt.figure(figsize=(6, 5))
    if power:
        scale = power_fit(data, power)
    else:
        [ scale, power ] = powerlaw_fit(data)
    exp_fit = exponential_fit(data)
    plot_data(data, (scale, power), exp_fit)
    # label axes
    plt.xlim(4, 18)
    plt.xlabel("Cell size (2.PI.R)", fontsize=14)
    #plt.ylim(45, 360)
    plt.ylabel(label, fontsize=14)
    plt.title(title, fontsize=18)
    plt.savefig('p-'+label+'.pdf', dpi=300)
    plt.close()


def plot_scaling(arg):
    L = [ 2*math.pi*d[1] for d in arg ]
    plot(zip(L, [ d[6] for d in arg ]), 'polymer', 2)
    plot(zip(L, [ d[5] for d in arg ]), 'microtubules', 2)
    plot(zip(L, [-d[7] for d in arg ]), 'ring-tension', [])
    plot(zip(L, [ d[8] for d in arg ]), 'radial-force', [])
    plot(zip(L, [ d[9] for d in arg ]), 'ring-length', 1)
    laplace(zip(L, [-d[7] for d in arg], [d[8] for d in arg]))


def clak(data, label):
    fig, ax = plt.subplots(figsize=(6, 5))
    # add mean:
    M = 0
    S = 0;
    for x, y in data:
        ax.plot(x, y, 'o', markersize=4, markerfacecolor='none', markeredgewidth=0.5, markeredgecolor=[0,0,1], alpha=0.7)
        M = M + y
        S = S + 1
    # plot mean:
    if S > 0:
        M = M / S
        L, R = ax.get_xlim()
        ax.plot([0, R], [M, M], '-', color='blue', linewidth=0.5, alpha=0.5)
    # label axes
    ax.set_xlabel("simulation", fontsize=14)
    #ax.set_ylim(45, 360)
    ax.set_ylabel(label, fontsize=14)
    ax.set_title(title, fontsize=18)
    plt.savefig('n-'+label+'.pdf', dpi=300)
    plt.close()


def plot_repeat(arg):
    N = [ int(d[0]) for d in arg ]
    clak(zip(N, [ d[6] for d in arg ]), 'polymer')
    clak(zip(N, [ d[5] for d in arg ]), 'microtubules')
    clak(zip(N, [-d[7] for d in arg ]), 'ring-tension')
    clak(zip(N, [ d[8] for d in arg ]), 'radial-force')
    clak(zip(N, [ d[9] for d in arg ]), 'ring-length')


#------------------------------------------------------------------------

def column_average(filename, start):
    """
    Calculate mean of values beyond line 'start'
    """
    inx = 0
    cnt = 0
    from operator import add
    with open(filename, 'r') as f:
        for line in f:
            strings = uncode(line).split()
            numbers = [ float(s) for s in strings ]
            inx = inx + 1
            if inx > start:
                cnt = cnt + 1
                if cnt > 1:
                    vec = list(map(add, vec, numbers))
                else:
                    vec = numbers
    return [ round(v/cnt, 3) for v in vec ]


def get_parameters(path):
    try:
        proc = subprocess.Popen(['grep', 'dimensions', path], stdout=subprocess.PIPE)
        data = uncode(proc.stdout.readline()).split()
        #print("config", data)
        res = [ float(data[2]), float(data[3]), float(data[4]) ]
        proc.stdout.close()
    except:
        res = [ 0, 0, 0 ]
    return res
    
    
#------------------------------------------------------------------------

def process(path):
    """ 
        This extracts parameters from the config file,
        and values from files 'ring.txt' or 'platelet.txt'
    """
    if path.startswith('run'):
        res = [ path[3:] ]
    else:
        res = [ path ]
    os.chdir(path)
    if not os.path.isfile('properties.cmo'):
        return []
    res.extend(get_parameters('config.cym'))    
    res.append('nan')
    # get data generated by "report time platelet verbose=0 > platelet.txt"
    if os.path.isfile('platelet.txt'):
        data = column_average('platelet.txt', 0)
        res.extend(data[1:])
    elif os.path.isfile('ring.txt'):
        res.extend(column_average('ring.txt', 0))
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
            if data:
                print(format_line(data))
                if nb_columns != len(data):
                    if nb_columns == 0:
                        nb_columns = len(data)
                    else:
                        sys.stderr.write("Error: data size mismatch in %s\n" % p)
                        break
                res.append(data)
        os.chdir(cdir)
        if data:
            plot_scaling(res)
            #plot_repeat(res)
        else:
            sys.stderr.write("  Error: no data!\n")
            sys.exit()


#------------------------------------------------------------------------

if __name__ == "__main__":
    if len(sys.argv) < 2 or sys.argv[1].endswith("help"):
        print(__doc__)
    else:
        main(sys.argv[1:])

