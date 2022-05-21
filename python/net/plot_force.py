#!/usr/bin/env python3
#
# Make a plot using matplotlib
#
# F. Nedelec, February 2016
#


"""
    Make plots
    
Syntax:
    
    plot_force.py DIRECTORY_PATH
    
Description:
    
"""

import sys, os, math, subprocess, random

try:
    import numpy as np
except ImportError:
    print("  Error: could not load numpy in python " + sys.version)
    sys.exit()

try:
    import matplotlib
except ImportError:
    print("  Error: could not load matplotlib in python " + sys.version)
    sys.exit()


import matplotlib.pyplot as plt
#matplotlib.use('SVG')

fts = 14


def uncode(arg):
    try:
        if isinstance(arg, unicode):
            return str(arg.decode('utf-8'))
    except:
        pass
    return arg


def plot_force(scale, data):
    #fig = plt.figure(figsize=(1, 1))
    #ax = fig.add_axes([0, 0, 1, 1])
    fig = plt.figure(figsize=(7, 5))
    xdata = np.array(scale)
    for k in data.keys():
        val = data[k]
        G = random.uniform(0,1)
        B = random.uniform(0,1)
        for v in val:
            R = random.uniform(0,1)
            ydata = np.array(v[0])
            plt.plot(xdata, ydata, '-', linewidth=0.5, color=(R, G, B))
    #print(k, len(ydata))
    #plt.xlim(0, 100)
    #plt.ylim(0, ymax)
    #plt.ylim(15, 30)
    #ax.set_xticks([])
    #ax.set_yticks([])
    plt.xlabel('Force (pN)', fontsize=fts)
    plt.ylabel('Count', fontsize=fts)
    plt.title('Force distribution', fontsize=fts)
    fig.tight_layout()
    #ax.grid(True)


def get_force(filename):
    """
    Work in current directory
    """
    if not os.path.isfile(filename):
        subprocess.call(['report2', 'couple:force'], stdout=open(filename, 'w'))
    scale = []
    data = { }
    with open(filename, 'r') as f:
        for line in f:
            spl = uncode(line).split()
            if len(spl) < 2:
                pass
            elif spl[0] == '%':
                pass
            elif len(spl) > 2:
                k = spl[0].strip()
                v = spl[1:-1]
                if k == 'scale':
                    scale = v
                else:
                    if k in data:
                        data[k].append([v])
                    else:
                        data[k] = [[v]]
    return scale, data


def parse(dirpath):
    """
        Work in current directory
        """
    scale, data = get_force('couple_force.txt')
    plot_force(scale, data)
    plt.savefig('couple_force', dpi=150)
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
        parse('.')
    else:
        cdir = os.getcwd()
        for p in paths:
            os.chdir(p)
            sys.stdout.write('- '*32+p+"\n")
            try:
                parse(p)
            except Exception as e:
                out.write("Error: %s\n" % repr(e))
            os.chdir(cdir)


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1].endswith("help"):
        print(__doc__)
    else:
        main(sys.argv[1:])

