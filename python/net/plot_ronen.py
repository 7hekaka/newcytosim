#!/usr/bin/env python3
#
# A script to plot for project with Ronen Zaidel-Bar
#
# F. Nedelec, Strasbourg, 16.12.2021, 8.1.2022


"""
    Plot the radius of the network as a function of time.
    This relies on 'reportN' to produce data in 'mom.txt'
    
Syntax:
    
    plot_ronen.py DIRECTORY_PATH
    
Description:
    
"""

#font size:
fts = 14
do_plot = 1
add_fit = 1
results = []

import sys, os, math, subprocess
try:
    import matplotlib
    #matplotlib.use('SVG')
    import matplotlib.pyplot as plt
except:
    do_plot = 0
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


def fit_curve(time, data):
    """
        Fit exponential and save results
    """
    (A, B) = exponential_fit(zip(time, data))
    return (A, -B)


def nice_plot(time, data):
    """
        Plot size as a function of time
    """
    fig = plt.figure(figsize=(5, 4))
    plt.plot(time, data, 'b-', linewidth=7)
    # add horizontal bar at starting size:
    h = data[0]
    plt.plot([min(time), max(time)], [h, h], 'k-', linewidth=1)
    #plt.xlim(0, 100)
    if add_fit:
        (A, B) = exponential_fit(zip(time, data))
        fit = [ A*math.exp(B*t) for t in time ]
        plt.plot(time, fit, 'w.', markersize=7)
    plt.ylim(0, 7)
    plt.xlabel('Time (s)', fontsize=fts)
    plt.ylabel('Radius (um)', fontsize=fts)
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
            M = float(s[8]) # sum of variances in X, Y, Z
        elif len(s) == 10 and s[1].isalpha():
            T = []
            M = []
            tim.append(float(s[0]))
            mom.append(float(s[9]))
    return tim, mom


def get_radius(file):
    """
        Extract radius from moment variance
    """
    T, M = get_moment(file)
    if not T:
        raise Exception("Could not find time information")
    R = [ math.sqrt(2*x) for x in M ]
    #S = [ 2*math.pi*x for x in M ]
    return T, R


def process(dirpath, filename):
    """
        Process given directory
    """
    os.chdir(dirpath)
    if not os.path.isfile(filename):
        args = ['reportN', 'fiber:moment']
        subprocess.call(args, stdout=open(filename, 'w'))
    with open(filename, 'r') as f:
        time, data = get_radius(f)
        if do_plot:
            nice_plot(time, data)
            plt.savefig('size.png', dpi=150)
            #plt.show()
            plt.close()
        if add_fit:
            (A, B) = fit_curve(time, data)
            results.append([dirpath, A, B])


def print_results(filename=''):
    if filename:
        file = open(filename, 'w')
    else:
        file = sys.stdout
    for i in results:
        S = round(i[1]*math.exp(-5*i[2]), 3)
        file.write("%s %f %f\n" % (i[0], S, i[2]))
    if filename:
        file.close()


def read_results(filename):
    """
        Read numeric data from file
    """
    data = []
    with open(filename, 'r') as f:
        for line in f:
            if len(line) > 5 and line[0]!='%':
                s = line.split()
                data.append([s[0], float(s[1]), float(s[2])])
    return data


def master_plot(data):
    """
        Read numeric data from file
    """
    P, A, B = zip(*data)
    B1 = B[0:len(B)-1:2]
    B2 = B[1:len(B):2]
    fig = plt.figure(figsize=(5, 4))
    plt.plot(B1, B2, 'bo', markersize=7)
    # add diagonal:
    h = math.ceil(10*max(B))*0.1
    plt.plot([0, h], [0, h], 'k-', linewidth=1)
    plt.xlim(0, h)
    plt.ylim(0, h)
    plt.xlabel('Contraction Rate (mod=0)', fontsize=fts)
    plt.ylabel('Contraction Rate (mod=1)', fontsize=fts)
    plt.title('Rate correlation', fontsize=fts)
    fig.tight_layout()
    plt.savefig('result.png', dpi=150)

#-----------------------------------------------------------------------------

def main(args):
    paths = []
    files = []
    global results
    for arg in args:
        if os.path.isdir(arg):
            paths.append(arg)
        elif os.path.isfile(arg):
            files.append(arg)
        else:
            sys.stderr.write("  Error: unexpected argument `%s'\n" % arg)
            sys.exit()
    
    if not paths:
        if files:
            results = read_results(files[0])
        else:
            process('.', 'mom.txt')
    else:
        cdir = os.getcwd()
        for p in paths:
            sys.stdout.write('- '*32+p+"\n")
            process(p, 'mom.txt')
            os.chdir(cdir)
    if do_plot:
        master_plot(results)
    print_results()


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1].endswith("help"):
        print(__doc__)
    else:
        main(sys.argv[1:])

