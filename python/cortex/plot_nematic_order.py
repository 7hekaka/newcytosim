#!/usr/bin/env python3
#
# Plot nematic order, calculated by `report fiber:nematic`
#
# F. Nedelec, Strasbourg, 06-07.12.2024


"""
Description:
    Plot data calculated by `report fiber:nematic verbose=0 > order.txt`

Syntax:
    plot_nematic_order.py DIRECTORY_PATH

"""

# time cut off (seconds) above which data is averaged
cutoff = 4000

#font size:
fts = 14
# file format
format = 'png'

import sys, os, math, subprocess
import matplotlib
import matplotlib.pyplot as plt
#matplotlib.use('SVG')


def uncode(arg):
    try:
        return arg.decode('utf-8')
    except:
        return arg


def retreive_data(file):
    """
        Retreive data from file
    """
    T = []
    D = []
    ts = 0
    for line in file:
        v = uncode(line).split()
        t = 0
        d = []
        if len(v) < 2:
            pass
        elif v[0] == '%' and v[1] == "time":
            t = float(v[2])
        elif len(v) == 7:
            t = float(v[1])
            d = float(v[3])
        if t < ts:
            sys.stderr.write(f'Warning: discarding {len(T)} duplicated datalines\n')
            T = []
            D = []
            ts = 0
            continue
        if d:
            T.append(t)
            D.append(d)
            ts = t
    return T, D


def plot_data(X, Y, M, name):
    """
        Plot data as a function of time
    """
    Xinf = math.floor(min(X)/100)*100
    Xsup = math.ceil(max(X)/100)*100
    Yinf = 0.2
    Ysup = math.ceil(max(Y))
    fig = plt.figure(figsize=(4, 3))
    plt.scatter(X, Y, label='nematic order', marker='o', s=8)
    plt.plot([Xinf, Xsup], [M, M], linewidth=1.0)
    plt.xlim(Xinf, Xsup)
    plt.ylim(Yinf, Ysup)
    plt.xlabel('Time (s)', fontsize=fts)
    plt.ylabel('Order', fontsize=fts)
    plt.title('Nematic order', fontsize=fts)
    #plt.legend()
    fig.tight_layout()
    if format == 'svg':
        plt.savefig('order.svg', format='svg', dpi=150)
    else:
        plt.savefig('order.png', dpi=150)
    #plt.show()
    plt.close()


def process(dirpath):
    """
        Generate a plot in given directory
    """
    filename = "order.txt"
    if not os.path.isfile(filename):
        if not os.path.isfile("objects.cmo"):
            sys.stderr.write(f'Warning: missing file {dirpath}/objects.cmo\n')
            return
        else:
            # attempt to generate report file:
            args = ["report3", "fiber:nematic", "verbose=0"]
            subprocess.call(args, stdout=open(filename, 'w'))
    res = 0
    with open(filename, 'r') as f:
        T, D = retreive_data(f)
        avg = math.nan
        if T:
            #print(T, D)
            # calculate mean length for data above cutoff:
            val = [ x for t,x in zip(T,D) if t > cutoff ]
            if val:
                avg = sum(val) / len(val)
            plot_data(T, D, avg, dirpath)
        else:
            sys.stderr.write(f'Warning: no data in {dirpath}/{filename}\n')
    print(f'{dirpath} {avg:6.3f}')


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
            sys.stderr.write(f"Error: unexpected argument `{arg}`\n")
            sys.exit()
    if not paths:
        process('.')
    else:
        cdir = os.getcwd()
        for p in paths:
            os.chdir(p)
            process(p)
            os.chdir(cdir)


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1].endswith("help"):
        print(__doc__)
    else:
        main(sys.argv[1:])

