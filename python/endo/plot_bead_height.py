#!/usr/bin/env python3
#
# Plot bead Z position, calculated by `report solid`
#
# F. Nedelec, Isola di Servolo, Venice 11.09.2023


"""
Description:
    Plot bead Z-position, calculated by:
    report3 solid verbose=0 > bead.txt
    
Syntax:
    plot_bead_heigth.py DIRECTORY_PATH

To finish the plot:

"""

#font size:
fts = 14

import sys, os, math, subprocess
import matplotlib
import matplotlib.pyplot as plt
#matplotlib.use('SVG')


def uncode(arg):
    try:
        return arg.decode('utf-8')
    except:
        return arg


def plot_data(X, Y, name, dots=0):
    """
        Plot fiber length as function of frame
    """
    fig = plt.figure(figsize=(4, 3))
    if dots:
        plt.plot(X, Y, linewidth=0, marker='.', markersize=dots)
        plt.xlabel('SIMUL', fontsize=fts)
    else:
        plt.plot(X, Y, linewidth=4)
        plt.xlabel('Time', fontsize=fts)
    plt.xlim(0, math.ceil(max(X)))
    plt.ylim(-0.030, 0.030) #math.ceil(max(Y)))
    plt.ylabel(r'Z-position ($\mu m$)', fontsize=fts)
    plt.title('Bead height', fontsize=fts)
    #plt.legend()
    fig.tight_layout()
    plt.savefig('bead.png', dpi=75*(1+dots))
    #plt.show()
    plt.close()


def get_data(file):
    """
        Retreive data from file
    """
    X = []
    Y = []
    for line in file:
        s = uncode(line).split()
        if len(s) < 2:
            pass
        elif s[0] in ('%', '#'):
            pass
        elif len(s) == 12:
            T = float(s[0])
            Z = float(s[5])
            #print(T,Z)
            X.append(T)
            Y.append(Z)
    return X, Y


def process(dirpath):
    """
        Process given directory
    """
    os.chdir(dirpath)
    filename = 'bead.txt'
    if not os.path.isfile(filename):
        args = ['report3', 'time', 'solid', 'verbose=0']
        subprocess.call(args, stdout=open(filename, 'w'))
    res = 0
    with open(filename, 'r') as f:
        X, Y = get_data(f)
        #print(D, N)
    plot_data(X, Y, dirpath)
    # record position at last time point:
    return X[-1], Y[-1]


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
        process('.')
    else:
        i = 0
        X = []
        Y = []
        cdir = os.getcwd()
        for p in paths:
            T, Z = process(p)
            i = i+1
            X.append(i)
            Y.append(Z)
            os.chdir(cdir)
        plot_data(X, Y, cdir, 7)

if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1].endswith("help"):
        print(__doc__)
    else:
        main(sys.argv[1:])

