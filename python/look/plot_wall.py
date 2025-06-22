#!/usr/bin/env python3
#
# Make a plot using matplotlib
#
# F. Nedelec, 13-15 November 2013, 22.06.2025
#
# How I installed matplotlib in Jan 2015 (Mac osx 10.10):
#
# brew install python3 --framework
# pip3 install nose
# pip3 install matplotlib
#
# Alternative from source:
# git clone git://github.com/matplotlib/matplotlib.git
# cd matplotlib
# python3 setup.py build
# sudo python3 setup.py install
# python3 setup.py clean
#


"""
    Make plots
    
Syntax:
    
    plot_stuff.py DIRECTORY-PATH
    
Description:
    
"""

import sys, os, math, subprocess
import matplotlib.pyplot as plt
import read_config

fts = 14
out = sys.stderr


def load_data(filename):
    """
    read line-by-line and convert to numerical value is possible
    """
    res = []
    with open(filename, 'r') as f:
        for line in f:
            if line and line[0] != '%':
                L = []
                for v in line.split():
                    try:
                        v = float(v)
                        if v == int(v):
                            v = int(v)
                    except:
                        pass
                    L.append(v)
                if L:
                    res.append(L)
    # reformat data in columns:
    return list(zip(*res))


def image(x, p):
	while x > p:
		x = x - 2*p
	while x < -p:
		x = x + 2*p
	return x


def plot_diameter(exp, sim):
    plt.figure(figsize=(3.84, 2.56), dpi=100)
    plt.plot(exp[0], exp[1], 'ko')
    plt.plot(sim[0], sim[1], 'b-', linewidth=4)
    plt.xlabel('Time (s)', fontsize=fts)
    plt.ylabel(r'Diameter ($\mu m$)', fontsize=fts)
    #plt.title('Closure', fontsize=fts)
    plt.tight_layout()
    #plt.xlim([-5, 5])
    #plt.grid(True)


def plot_positions(xdata, xlim, xinfo):
    plt.figure(figsize=(3.84, 2.56), dpi=100)
    n, bins, patches = plt.hist(xdata, 20)
    ax = plt.gca()
    ax.yaxis.set_visible(False)
    plt.xlabel('Position (%s)' % xinfo, fontsize=fts)
    plt.ylabel('Amount', fontsize=fts)
    plt.title('Histogram', fontsize=fts)
    plt.tight_layout()
    plt.xlim(xlim)
    #plt.grid(True)


def parse(dirpath):
    """
    Extract data in current directory
    """
    pile = read_config.parse('config.cym')
    shape = read_config.get_value(pile, ['set', 'space', 'cell', 'shape'])
    print("shape is [",shape,"]")
    if shape.startswith('wall'):
        if not os.path.isfile('radius.txt'):
            subprocess.call(['reportW', 'space'], stdout=open('radius.txt', 'w'), stderr=None)
        data = load_data('radius.txt')
        if len(data) > 0:
            tim = data[3]
            rad = [ 2*x for x in data[4] ]
            exp = load_data('/Users/nedelec/tmp/closure_diameter.txt')
            plot_diameter(exp, [tim, rad])
            plt.savefig('diameter.png')
            plt.close()
        #
        if not os.path.isfile('beads.txt'):
            subprocess.call(['reportW', 'bead', 'frame=1999'], stdout=open('beads.txt', 'w'), stderr=None)
        data = load_data('beads.txt')
        data = [ math.atan2(data[3][n], data[2][n]) for n in range(len(data[0])) ]
        if len(data) > 0:
            plot_positions(data, [-3.1416, 3.1416], 'radian')
            plt.savefig('position.png')
            plt.close()

    if shape.startswith('strip'):
        if not os.path.isfile('beads.txt'):
            subprocess.call(['reportW', 'bead', 'frame=999'], stdout=open('beads.txt', 'w'), stderr=None)
        data = load_data('beads.txt')
        if len(data) > 0:
            data = [ image(x, 5) for x in data[2] ]
            plot_positions(data, [-5, 5], 'um')
            plt.savefig('position.png')
            plt.close()


#------------------------------------------------------------------------

def main(args):
    paths = []
    
    for arg in args:
        if os.path.isdir(arg):
            paths.append(arg)
        else:
            out.write("  Error: unexpected argument `%s'\n" % arg)
            sys.exit()
    
    if not paths:
        parse('.')
    else:
        cdir = os.getcwd()
        for p in paths:
            os.chdir(p)
            sys.stdout.write('- '*32+'\n')
            parse(p)
            os.chdir(cdir)

#------------------------------------------------------------------------

if __name__ == "__main__":
    if len(sys.argv) < 2 or sys.argv[1].endswith("help"):
        print(__doc__)
    else:
        main(sys.argv[1:])

