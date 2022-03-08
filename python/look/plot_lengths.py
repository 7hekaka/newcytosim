#!/usr/bin/env python3
#
# Make a plot using matplotlib
#
# F. Nedelec, Strasbourg, 08.03.2022


"""
Description:
    Plot the histogram of filament length distribution
    This relies on 'reportN' to produce data in 'len.txt'
    
Syntax:
    plot_lengths.py DIRECTORY_PATH

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


def plot_histograms(time, data, scale, name):
    """
        Plot surface as a function of time
    """
    fig = plt.figure(figsize=(4, 3.5))
    Y = data[-1]
    W = (scale[1]-scale[0])*0.8
    plt.bar(scale, Y, W, label=name)
    sup = 10 * math.ceil((min(Y)+max(Y))/10)
    plt.ylim(0, sup)
    plt.xlabel('Length (um)', fontsize=fts)
    plt.ylabel('Count', fontsize=fts)
    plt.title('Length histogram', fontsize=fts)
    plt.legend()
    fig.tight_layout()


def get_lengths(file):
    """
        Retreive length histrogram from file
    """
    T = []
    D = []
    S = []
    F = 'fiber'
    for line in file:
        s = uncode(line).split()
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
    if 1:
        filename='len.txt'
        if not os.path.isfile(filename):
            subprocess.call(['reportN', 'fiber:histogram', 'interval=0.5,10'], stdout=open(filename, 'w'))
        with open(filename, 'r') as f:
            T, D, S, F = get_lengths(f)
            plot_histograms(T, D, S, F)
    plt.savefig('lengths.png', dpi=150)
    #plt.show()
    plt.close()

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

