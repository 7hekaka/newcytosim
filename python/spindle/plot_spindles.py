#!/usr/bin/env python3
#
# Masterplot for Arabidopsis Spindle project
#
# F. Nedelec, Cambridge 23.06.2023


"""
Description:
    Make different master plots
    
Syntax:
    plot_spindles.py

"""

#font size:
fts = 14

import sys, os, math

import matplotlib
import matplotlib.pyplot as plt
#matplotlib.use('SVG')


def uncode(arg):
    try:
        return arg.decode('utf-8')
    except:
        return arg


def read_data_file(path):
    """
        get data from file in column format
    """
    res = []
    with open(path, 'r') as file:
        for line in file:
            code = uncode(line).split()
            data = []
            #print(code)
            if code[0].startswith('run'):
                data = [code[0]]
            else:
                data = [code[0]]
            for i in code[1:]:
                data.append(float(i))
            res.append(data)
    res = list(zip(*res))
    #print(res)
    return res


def plot_spindle_length(X, L):
    """
        Plot surface as a function of time
    """
    fig = plt.figure(figsize=(4, 3))
    plt.plot(X, L, marker='o', markersize=4, linewidth=0, markeredgecolor='none')
    plt.xlim(1, math.ceil(max(X)))
    plt.ylim(math.floor(min(L)), math.ceil(max(L)))
    plt.xlabel('Augmin level', fontsize=fts)
    plt.ylabel('Pole-to-pole distance (um)', fontsize=fts)
    plt.title('Spindle Length', fontsize=fts)
    plt.legend()
    fig.tight_layout()
    plt.savefig('spindle_length.png', dpi=150)



def plot_fiber_lengths(X, L0, L1, L2):
    """
        Plot surface as a function of time
    """
    fig = plt.figure(figsize=(4, 3))
    plt.plot(X, L0, marker='o', markersize=4, linewidth=0, color='green')
    plt.plot(X, L1, marker='o', markersize=4, linewidth=0, color='orange')
    plt.plot(X, L2, marker='o', markersize=4, linewidth=0, color='blue')
    plt.xlim(1, math.ceil(max(X)))
    mL = min(min(L0), min(L1), min(L2))
    xL = max(max(L0), max(L1), max(L2))
    plt.ylim(math.floor(mL), math.ceil(xL))
    plt.xlabel('Augmin level', fontsize=fts)
    plt.ylabel('Fiber length (um)', fontsize=fts)
    plt.title('Mean fiber Lengths', fontsize=fts)
    plt.legend()
    fig.tight_layout()
    plt.savefig('fiber_length.png', dpi=150)


#------------------------------------------------------------------------

def main(args):
    pam = read_data_file('parameters.txt')
    len = read_data_file('spindle_length.txt')
    mtl = read_data_file('fiber_length.txt')
    plot_spindle_length(pam[1], len[1])
    plot_fiber_lengths(pam[1], mtl[1], mtl[2], mtl[3])


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1].endswith("help"):
        print(__doc__)
    else:
        main(sys.argv[1:])

