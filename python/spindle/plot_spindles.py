#!/usr/bin/env python3
#
# Masterplot for Arabidopsis Spindle project
#
# F. Nedelec, Cambridge 23.06.2023


"""
Description:

    Make different master plots for spindle simulation analysis

Usage:
    get_parameters.py run* > parameters.txt
    plot_spindle_length.py run* > spindle_length.txt
    plot_fiber_length.py run* > fiber_length.txt

    plot_spindles.py

23.06.2023
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
        read data from file and return in column format
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
    return res


def plot_spindle_length(X, Y):
    """
        Plot surface as a function of time
    """
    fig = plt.figure(figsize=(4, 3))
    plt.plot(X, Y, marker='o', markersize=4, linewidth=0, markeredgecolor='none')
    plt.xlim(1, math.ceil(max(X)))
    plt.ylim(math.floor(min(Y)), math.ceil(max(Y)))
    plt.xlabel('Augmin level', fontsize=fts)
    plt.ylabel(r'Pole-to-pole distance ($\mu m$)', fontsize=fts)
    plt.title('Spindle Length', fontsize=fts)
    plt.legend()
    fig.tight_layout()
    plt.savefig('spindle_length.png', dpi=150)



def plot_fiber_lengths(X, Y0, Y1, Y2):
    """
        Plot surface as a function of time
    """
    fig = plt.figure(figsize=(4, 3))
    plt.plot(X, Y0, marker='o', markersize=4, linewidth=0, color='green')
    plt.plot(X, Y1, marker='o', markersize=4, linewidth=0, color='orange')
    plt.plot(X, Y2, marker='o', markersize=4, linewidth=0, color='blue')
    plt.xlim(1, math.ceil(max(X)))
    mY = min(min(Y0), min(Y1), min(Y2))
    xY = max(max(Y0), max(Y1), max(Y2))
    plt.ylim(math.floor(mY), math.ceil(xY))
    plt.xlabel('Augmin level', fontsize=fts)
    plt.ylabel(r'Fiber length ($\mu m$)', fontsize=fts)
    plt.title('Mean fiber Lengths', fontsize=fts)
    plt.legend()
    fig.tight_layout()
    plt.savefig('fiber_length.png', dpi=150)


def plot_fiber_count(X, Y0, Y1, Y2):
    """
        Plot surface as a function of time
    """
    fig = plt.figure(figsize=(4, 3))
    plt.plot(X, Y0, marker='o', markersize=4, linewidth=0, color='green')
    plt.plot(X, Y1, marker='o', markersize=4, linewidth=0, color='orange')
    plt.plot(X, Y2, marker='o', markersize=4, linewidth=0, color='blue')
    plt.xlim(1, math.ceil(max(X)))
    mY = min(min(Y0), min(Y1), min(Y2))
    xY = max(max(Y0), max(Y1), max(Y2))
    plt.ylim(math.floor(mY), math.ceil(xY))
    plt.xlabel('Augmin level', fontsize=fts)
    plt.ylabel('Fiber count', fontsize=fts)
    plt.title('Mean fiber counts', fontsize=fts)
    plt.legend()
    fig.tight_layout()
    plt.savefig('fiber_count.png', dpi=150)


#------------------------------------------------------------------------

def main(args):
    try:
        pam = read_data_file('parameters.txt')
    except FileNotFoundError:
        print("please run `get_parameters.py run* > parameters.txt'");
        return
    try:
        spi = read_data_file('spindle_length.txt')
    except FileNotFoundError:
        print("please run `plot_spindle_length.py run* > spindle_length.txt'");
        return
    try:
        fib = read_data_file('fiber_length.txt')
    except FileNotFoundError:
        print("please run `plot_fiber_length.py run* > fiber_length.txt'");
        return
    plot_spindle_length(pam[1], spi[1])
    plot_fiber_count(pam[1], fib[1], fib[3], fib[5])
    plot_fiber_lengths(pam[1], fib[2], fib[4], fib[6])


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1].endswith("help"):
        print(__doc__)
    else:
        main(sys.argv[1:])

