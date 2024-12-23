#!/usr/bin/env python3
#
# Make a plot using matplotlib
#
# F. Nedelec, Strasbourg, 08.03.2022, Cambridge 15.06.2023


"""
Description:
    Plot the histogram of filament length distribution
    This relies on 'report3' to produce data in 'len.txt'
    
Syntax:
    plot_lengths.py DIRECTORY_PATH

"""

# size of font used in plots:
fts = 14
# file format: `png` or `svg`
format = 'png'

import sys, os, math, subprocess
import matplotlib
import matplotlib.pyplot as plt


def split_line(arg):
    """
        Not sure if this is useful... was needed in older python
    """
    try:
        return arg.decode('utf-8').split()
    except:
        return arg.split()


def plot_histogram(time, data, scale, name):
    """
        Plot surface as a function of time
    """
    fig = plt.figure(figsize=(4, 3))
    Y = data[-1]
    W = (scale[1]-scale[0])*0.8
    plt.bar(scale, Y, W, label=name)
    plt.xlim(0, max(scale))
    sup = 10 * math.ceil((min(Y)+max(Y))/10)
    plt.ylim(0, sup)
    plt.xlabel(r'Length ($\mu m$)', fontsize=fts)
    plt.ylabel('Count', fontsize=fts)
    plt.title('Length histogram', fontsize=fts)
    plt.legend()
    fig.tight_layout()


def get_lengths(file):
    """
        Retreive length histogram from file
    """
    T = []
    D = []
    S = []
    F = 'fiber'
    for line in file:
        s = split_line(line)
        if len(s) < 2:
            pass
        elif s[0] == '%':
            if s[1] == "time":
                T.append(float(s[2]))
        elif s[0] == "scale":
            S = [ float(x) for x in s[1:] ]
        else:
            F = s[0]
            H = [ float(x) for x in s[1:] ]
            D.append(H)
    return T, D, S, F


def process(dirpath):
    """
        Process given directory
    """
    os.chdir(dirpath)
    filename='len.txt'
    args = ['report3', 'fiber:histogram', 'interval=0.1,32', 'verbose=0']
    subprocess.call(args, stdout=open(filename, 'w'))
    with open(filename, 'r') as f:
        T, D, S, F = get_lengths(f)
        plot_histogram(T, D, S, F)
    #plt.show()
    if format == 'svg':
        plt.savefig('lengths.svg', format='svg', dpi=150)
    else:
        plt.savefig('lengths.png', dpi=150)
    plt.close()

#------------------------------------------------------------------------

def main(args):
    global format
    paths = []
    
    for arg in args:
        if os.path.isdir(arg):
            paths.append(arg)
        elif arg in ('png', 'svg' ):
            format = arg
        else:
            sys.stderr.write("  Error: unexpected argument `%s'\n" % arg)
            sys.exit()
    
    if not paths:
        process('.')
    else:
        cdir = os.getcwd()
        for p in paths:
            sys.stdout.write('- '*32+p+"\n")
            process(p)
            os.chdir(cdir)


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1].endswith("help"):
        print(__doc__)
    else:
        main(sys.argv[1:])

