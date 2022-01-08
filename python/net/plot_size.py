#!/usr/bin/env python3
#
# Make a plot using matplotlib
#
# F. Nedelec, February 2016
# Strasbourg, 16 December 2021


"""
    Plot the radius of the network as a function of time.
    This relies on 'reportN' to produce data in 'mom.txt'
    
Syntax:
    
    plot_size.py DIRECTORY_PATH
    
Description:
    
"""

#font size:
fts = 14
add_fit = 1

import sys, os, math, subprocess
import matplotlib
import matplotlib.pyplot as plt
#matplotlib.use('SVG')
try:
    from pyned import exponential_fit
except:
    add_fit = 0


def uncode(arg):
    try:
        return arg.decode('utf-8')
    except:
        return arg


def prune_values(time, data):
    """
        Clean dataset by removing infinite values
    """
    i = len(data)-1
    while i >= 0:
        if not math.isfinite(data[i]):
            data.pop(i)
            time.pop(i)
        i-=1


def nice_plot(time, data):
    """
        Plot surface as a function of time
    """
    fig = plt.figure(figsize=(5, 4))
    plt.plot(time, data, 'b-', linewidth=7)
    # add horizontal bar at starting size:
    h = data[0]
    plt.plot([min(time), max(time)], [h, h], 'k-', linewidth=1)
    # add exponential fit to data:
    if add_fit:
        (A, B) = exponential_fit(zip(time, data))
        fit = [ A*math.exp(B*t) for t in time ]
        plt.plot(time, fit, 'w.', markersize=7)
    #plt.xlim(0, 100)
    plt.ylim(0, min(data)+max(data))
    plt.xlabel('Time (s)', fontsize=fts)
    plt.ylabel('Surface (um^2)', fontsize=fts)
    plt.title('Network size', fontsize=fts)
    fig.tight_layout()


def mini_plot(time, size):
    """
        Make small plot of size as a function of time
    """
    fig = plt.figure(figsize=(2, 2))
    ax = fig.add_axes([0, 0, 1, 1])
    plt.plot(time, data, 'b-', linewidth=3)
    h = data[0]
    plt.plot([min(time), max(time)], [h, h], 'k-', linewidth=1)
    #plt.xlim(0, 100)
    plt.ylim(0, min(data)+max(data))


def get_moment(file):
    """
        Get network size as a function of time
    """
    tim = []
    mom = []
    T = 0
    M = 0
    for line in file:
        s = uncode(line).split()
        if len(s) < 2:
            pass
        elif s[0] == '%':
            if s[1] == "start" or s[1] == "time":
                T = float(s[2])
            elif s[1] == "moment":
                M = float(s[-1])
            elif s[1] == "end":
                tim.append(T)
                mom.append(M)
        elif len(s) == 9 and s[0].isalpha():
            M = float(s[8])
        elif len(s) == 10 and s[1].isalpha():
            T = []
            M = []
            tim.append(float(s[0]))
            mom.append(float(s[9]))
    return tim, mom


def get_size(file):
    """
        Extract surface from moment variance
    """
    T, M = get_moment(file)
    if not T:
        raise Exception("Could not find time information")
    S = [ 2*math.pi*x for x in M ]
    return T, S


def process(dirpath):
    """
        Process given directory
    """
    os.chdir(dirpath)
    if 1:
        filename='mom.txt'
        if not os.path.isfile(filename):
            subprocess.call(['reportN', 'fiber:moment'], stdout=open(filename, 'w'))
        with open(filename, 'r') as f:
            time, data = get_size(f)
            nice_plot(time, data)
    if 0:
        filename = 'cross.txt'
        if not os.path.isfile(filename):
            subprocess.call(['reportN', 'fiber:intersection'], stdout=open(filename, 'w'))
        with open(filename, 'r') as f:
            time, data = get_size(f)
            prune_values(time, data)
            plt.plot(time, data, 'k-', linewidth=2)
    if 0:
        filename = 'connectors.txt'
        if not os.path.isfile(filename):
            subprocess.call(['reportN', 'fiber:connector'], stdout=open(filename, 'w'))
        with open(filename, 'r') as f:
            time, data = get_size(f)
            prune_values(time, data)
            plt.plot(time, data, 'r-', linewidth=2)
    plt.savefig('size.png', dpi=150)
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

