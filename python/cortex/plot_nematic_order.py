#!/usr/bin/env python3
#
# Plot nematic order, calculated by `report fiber:nematic`
#
# F. Nedelec, Strasbourg, Chitry, 06-24.12.2024


"""
Description:
    Plot data calculated by `report fiber:nematic verbose=0 > order.txt`

Syntax:
    plot_nematic_order.py DIRECTORY_PATH

"""

# time cut off (seconds) above which data is averaged
cutoff = 4000

# size of font used in plots:
fts = 14
# file format: `png` or `svg`
format = 'png'
# name of data & figure files to be made
output = 'order'

import sys, os, math, subprocess
import matplotlib
import matplotlib.pyplot as plt


def uncode(arg):
    """
        Not sure if this is useful... was needed in older python
    """
    try:
        return arg.decode('utf-8')
    except:
        return arg


def retreive_data(filename):
    """
        Retreive data formatted as two columns, from file
    """
    T = []
    D = []
    ts = 0
    with open(filename, 'r') as file:
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


def plot_data(X, Y, M = math.nan):
    """
        Plot data as a function of time, if M is specified, add horizontal line
    """
    Xinf = math.floor(min(X)/100)*100
    Xsup = math.ceil(max(X)/100)*100
    Yinf = 0.2
    Ysup = math.ceil(max(Y))
    fig = plt.figure(figsize=(8, 6))
    plt.scatter(X, Y, label='Nematic Order', marker='o', s=3)
    if not math.isnan(M):
        plt.plot([Xinf, Xsup], [M, M], linewidth=1.0)
    plt.xlim(Xinf, Xsup)
    plt.ylim(Yinf, Ysup)
    plt.xlabel('Time (s)', fontsize=fts)
    plt.ylabel('Order', fontsize=fts)
    plt.title('Nematic Order', fontsize=fts)
    fig.tight_layout()
    #plt.legend()
    #plt.show()


def save_plot(name):
    if format == 'svg':
        plt.savefig(name+'.svg', format='svg')
    else:
        plt.savefig(name+'.png', dpi=150)
    plt.close()


def process(dirpath):
    """
        Generate a plot in given directory
    """
    filename = output+'.txt'
    if not os.path.isfile(filename):
        if not os.path.isfile("objects.cmo"):
            sys.stderr.write(f'Warning: missing file {dirpath}/objects.cmo\n')
            return
        else:
            # attempt to generate report file:
            print(f"plot_nematic_order.py makes {dirpath}/order.txt");
            args = ["report3", "fiber:nematic", "verbose=0"]
            subprocess.call(args, stdout=open(filename, 'w'))
    res = 0
    avg = math.nan
    T, D = retreive_data(filename)
    if T:
        #print(T, D)
        # calculate mean length for data above cutoff:
        val = [ x for t,x in zip(T,D) if t > cutoff ]
        if val:
            avg = sum(val) / len(val)
        plot_data(T, D, avg)
        save_plot(output);
    print(f'{dirpath} {avg:6.3f}')


def multi_plot(files):
    """
        Generate a plot in given directory
    """
    f = files[0]
    T, D = retreive_data(f)
    plot_data(T, D)
    for f in files[1:]:
        T, D = retreive_data(f)
        if T:
            plt.scatter(T, D, label='Nematic Order', marker='o', s=3)
    save_plot(output);


#------------------------------------------------------------------------

def main(args):
    global format, output
    paths = []
    files = []
    for arg in args:
        if os.path.isdir(arg):
            paths.append(arg)
        elif arg.endswith('.png') or arg.endswith('.svg'):
            [output, _, format] = arg.partition('.')
        elif os.path.isfile(arg):
            files.append(arg)
        elif arg in ('png', 'svg' ):
            format = arg
        else:
            sys.stderr.write(f"Error: unexpected argument `{arg}`\n")
            sys.exit()
    if files:
        multi_plot(files)
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

