#!/usr/bin/env python3
#
# Plot spindle length, calculated by `report spindle:length`
#
# F. Nedelec, Strasbourg, Cambridge 20--21.06.2023


"""
Description:
    Plot spindle length, calculated by:
    report time spindle:length > spindle_length.txt

Syntax:
    plot_spindle_length.py DIRECTORY_PATH

To get the augmin value:
    grep -H preconfig */config.cym > A
    sed "s/config.cym:%preconfig.augmin_source=/ /" A > augmin.txt
    cut -c 4-7 -c 9- augmin.txt > A
    cut -c 4- spindle_length.txt > L
    paste A L > AL.txt

"""

# time cut off
cutoff = 800

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
        elif len(s) == 10:
            T.append(float(s[0]))
            D.append(float(s[9]))
        elif len(s) > 7:
            T.append(ts)
            D.append(float(s[-1]))
    return T, D


def plot_data(X, Y, name):
    """
        Plot data as a function of time
    """
    fig = plt.figure(figsize=(4, 3))
    plt.plot(X, Y, label="pole-to-pole", linewidth=4.0)
    plt.xlim(0, math.ceil(max(X)/100)*100)
    plt.ylim(0, math.ceil(max(Y)))
    plt.xlabel('Time (s)', fontsize=fts)
    plt.ylabel(r'Length ($\mu m$)', fontsize=fts)
    plt.title('Spindle Length', fontsize=fts)
    plt.legend()
    fig.tight_layout()
    plt.savefig('spindle_length.png', dpi=150)
    #plt.show()
    plt.close()


def process(dirpath):
    """
        Process given directory
    """
    os.chdir(dirpath)
    filename = 'spindle_length.txt'
    if not os.path.isfile(filename):
        args = ['report3', 'time', 'spindle:length']
        subprocess.call(args, stdout=open(filename, 'w'))
    res = 0
    with open(filename, 'r') as f:
        T, L = get_data(f)
        #print(T, L)
        # calculate mean length for data above cutoff:
        LL = [ x for t,x in zip(T,L) if t > cutoff ]
        if LL:
            res = sum(LL) / len(LL)
        plot_data(T, L, dirpath)
    print(f'{dirpath} {res:6.3f}')


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

