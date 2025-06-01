#!/usr/bin/env python3
#
# Plot nematic order, calculated by `report fiber:nematic`
#
# F. Nedelec, Strasbourg, Chitry, 06-30.12.2024, 6.02.2025 30.05.2025


"""
Description:
    Plot data calculated by `report fiber:nematic verbose=0 > order.txt`

Syntax:
    plot_nematic_order.py DIRECTORY_PATH

"""

# time (seconds) above which data is averaged
threshold = 4000
# cut off for display of data on graphs
cutoff = 0

# size of marker in plots:
mks = 3
# size of font used in plots:
fts = 14
# file format: `png` or `svg`
format = 'png'
# name of data & figure files to be made
output = 'order'

import sys, os, math, subprocess, random
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
            v = line.split()
            t = ts
            d = []
            if len(v) < 2:
                pass
            elif v[0] == '%' and v[1] == "time":
                t = float(v[2])
            elif len(v) == 7:
                t = float(v[1])
                d = float(v[3])
            if t < ts:
                sys.stderr.write(f'Warning: discarding {len(T)} duplicated datalines in {filename}\n')
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
    Xinf = 0#math.floor(min(X)/100)*100
    Xsup = 8000#math.ceil(max(X)/100)*100
    Yinf = 0.2
    Ysup = math.ceil(max(Y))
    fig, axes = plt.subplots(figsize=(6, 4.5))
    axes.plot(X, Y, label='Nematic Order', marker='o', markersize=mks, color=(0.5,0,1))
    if not math.isnan(M):
        axes.plot([Xinf, Xsup], [M, M], linewidth=1.0)
    plt.xlim(Xinf, Xsup)
    plt.ylim(Yinf, Ysup)
    plt.xticks(fontsize=fts)
    plt.yticks([0.2, 0.4, 0.6, 0.8, 1.0], fontsize=fts)
    #plt.xlabel('Time (s)', fontsize=fts)
    #plt.ylabel('Order', fontsize=fts)
    #plt.title('Nematic Order', fontsize=fts)
    fig.tight_layout()
    #plt.legend()
    #plt.show()
    return fig, axes


def save_plot(name):
    if format == 'svg':
        plt.savefig(name+'.svg', format='svg')
    else:
        plt.savefig(name+'.png', dpi=150)
    plt.close()


def process(dirpath):
    """
        Generate data file and make plot in given directory
    """
    filename = output+'.txt'
    if not os.path.isfile(filename):
        if not os.path.isfile("objects.cmo"):
            sys.stderr.write(f'Warning: missing file {dirpath}/objects.cmo\n')
            return ([], [])
        else:
            # attempt to generate report file:
            sys.stderr.write(f"plot_nematic_order.py makes {dirpath}/order.txt");
            args = ["report3", "fiber:nematic", "verbose=0"]
            subprocess.call(args, stdout=open(filename, 'w'))
    print(f'{dirpath}', end=' ')
    T, D = retreive_data(filename)
    if T and D:
        #print(T, D)
        # calculate mean length for data above threshold:
        val = [ x for t,x in zip(T,D) if t > threshold ]
        if val:
            avg = sum(val) / len(val)
        else:
            avg = math.nan
        plot_data(T, D, avg)
        save_plot(output);
        print(f' {avg:6.3f} {max(T):9.0f}')
        return (T, D)
    else:
        return ([], [])


def redish_color():
    R = 0.5*random.random()
    G = 0.5*random.random()
    B = 0.5*random.random()
    return ( 1-R, G, 1-B )
    
def greenish_color():
    R = 0.5*random.random()
    G = 0.5*random.random()
    B = 0.5*random.random()
    return ( R, 1-G, B )
    
def blueish_color():
    R = 0.5*random.random()
    G = 0.5*random.random()
    B = 0.5*random.random()
    return ( R, G, 1-B )

#------------------------------------------------------------------------

def main(args):
    global format, output, cutoff
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
            [key, _, val] = arg.partition('=')
            if key == 'cutoff':
                cutoff = float(val)
            else:
                sys.stderr.write(f"Error: unexpected argument `{arg}`\n")
                sys.exit()
    if not paths and not files:
        paths = ['.']
    # generate/collect all data:
    data = []
    cdir = os.getcwd()
    for p in paths:
        os.chdir(p)
        T, D = process(p)
        if T and D:
            data.append((T,D))
        os.chdir(cdir)
    for f in files:
        T, D = retreive_data(f)
        if T and D:
            data.append((T,D))
    # make conbined plot:
    if data:
        T, D = data[0]
        Tx = [ x-cutoff for x in T ];
        fig, axes = plot_data(Tx, D)
        for T, D in data[1:5]:
            Tx = [ x-cutoff for x in T ];
            axes.plot(Tx, D, marker='o', markersize=mks, color=blueish_color())
        for T, D in data[5:]:
            Tx = [ x-cutoff for x in T ];
            axes.plot(Tx, D, marker='o', markersize=mks, color=greenish_color())
        save_plot(output);


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1].endswith("help"):
        print(__doc__)
    else:
        main(sys.argv[1:])

