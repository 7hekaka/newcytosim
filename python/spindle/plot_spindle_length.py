#!/usr/bin/env python3
#
# Plot spindle length, calculated by `report spindle:length`
#
# F. Nedelec, Strasbourg, Cambridge 20.06.2023


"""
Description:
    Plot spindle length, calculated by `report spindle:length`
    
Syntax:
    plot_spindle_length.py DIRECTORY_PATH

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


def plot_data(T, D, name):
    """
        Plot surface as a function of time
    """
    fig = plt.figure(figsize=(4, 3))
    plt.plot(T, D, label="pole-to-pole")
    plt.xlim(0, math.ceil(max(T)/100)*100)
    plt.ylim(0, math.ceil(max(D)))
    plt.xlabel('Time (s)', fontsize=fts)
    plt.ylabel('Length (um)', fontsize=fts)
    plt.title('Spindle Length', fontsize=fts)
    plt.legend()
    fig.tight_layout()


def get_data(file):
    """
        Retreive length from file
    """
    T = []
    D = []
    ts = 0
    for line in file:
        s = uncode(line).split()
        if len(s) < 2:
            pass
        elif s[0] == '%':
            if s[1] == "time":
                ts = float(s[2])
        elif len(s) > 7:
            T.append(ts)
            D.append(float(s[-1]))
    return T, D


def process(dirpath):
    """
        Process given directory
    """
    os.chdir(dirpath)
    filename = 'spindle_length.txt'
    args = ['report3', 'spindle:length', 'verbose=1']
    subprocess.call(args, stdout=open(filename, 'w'))
    res = 0
    with open(filename, 'r') as f:
        T, L = get_data(f)
        #print(T, L)
        plot_data(T, L, dirpath)
        # calculate mean length for data above 1000s:
        LL = [ x for t,x in zip(T,L) if t > 1000 ]
        res = sum(LL) / len(LL)
    plt.savefig('spindle_length.png', dpi=150)
    #plt.show()
    plt.close()
    print(f'{dirpath} {res}\n')


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
            process(p)
            os.chdir(cdir)


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1].endswith("help"):
        print(__doc__)
    else:
        main(sys.argv[1:])

