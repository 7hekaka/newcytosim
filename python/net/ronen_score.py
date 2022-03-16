#!/usr/bin/env python3
#
# A script to plot for project with Ronen Zaidel-Bar
#
# F. Nedelec, Strasbourg, 16.12.2021, 8-13.1.2022, 24.2.2022


"""
    Plot the radius of the network as a function of time.
    This relies on 'reportN' to produce data in 'mom.txt'
    All data sent to 'rates.txt'

Syntax:
    
    ronen_score.py DIRECTORY_PATHS
    
Description:
    
"""

#font size:
fts = 14
do_plot = 1
add_fit = 1
length = 0

# earliest time point
earliest = 0

import sys, os, math, subprocess
try:
    import matplotlib
    #matplotlib.use('SVG')
    import matplotlib.pyplot as plt
except:
    do_plot = 0
try:
    from pyned import exponential_fit, exponential_model_fit
except:
    add_fit = 0

#-------------------------------------------------------------------------------

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


def fit_curve(time, values):
    """
        Fit exponential and save results
    """
    data = list(zip(time, values))
    (A, B, C) = exponential_model_fit(data)
    (A, B) = exponential_fit(data)
    return (A, -B)


def nice_plot(time, data):
    """
        Plot size as a function of time
    """
    fig = plt.figure(figsize=(3.84, 3.84))
    plt.plot(time, data, 'b-', linewidth=7)
    # add horizontal bar at starting size:
    h = data[0]
    plt.plot([min(time), max(time)], [h, h], 'k-', linewidth=1)
    #plt.xlim(0, 100)
    if add_fit:
        (A, B) = exponential_fit(zip(time, data))
        fit = [ A*math.exp(B*t) for t in time ]
        txt = str(round(A,3)) + " exp("+str(round(B,5))+" t )"
        plt.plot(time, fit, 'w.', markersize=7, label=txt)
    plt.ylim(0, 10)
    plt.xlabel('Time (s)', fontsize=fts)
    plt.ylabel('Radius (um)', fontsize=fts)
    plt.legend(loc='upper right', fontsize=7)
    plt.title('Network size', fontsize=fts)
    fig.tight_layout()


def mini_plot(time, size):
    """
        Make small plot of size as a function of time
    """
    fig = plt.figure(figsize=(2.56, 1.92))
    ax = fig.add_axes([0, 0, 1, 1])
    plt.plot(time, data, 'b-', linewidth=3)
    h = data[0]
    plt.plot([min(time), max(time)], [h, h], 'k-', linewidth=1)
    #plt.xlim(0, 100)
    plt.ylim(0, min(data)+max(data))


def get_moment(file):
    """
        Get fiber's position variance as a function of time
    """
    res = []
    T = []
    M = []
    L = []
    for line in file:
        s = uncode(line).split()
        if len(s) < 2:
            pass
        elif s[0] == '%':
            if s[1] == "time":
                if T and M:
                    res.append([T, M, L])
                    M = []
                    L = []
                T = float(s[2])
            elif s[1] == "end":
                if T and M:
                    res.append([T, M, L])
                T = []
                M = []
                L = []
        elif len(s) == 9 and s[0].isalpha():
            M = float(s[8]) # sum of variances in X, Y, Z
            L = float(s[1])
    return zip(*res)


def print_data(data):
    for t, m, l in data:
        print("%.2f : %.2f %.2f" % (t, m, l))


def get_radius(file):
    """
        Compute radius of network from moment variance
    """
    global earliest, length
    T, M, L = get_moment(file)
    length = max(L) #sum(L) / len(L)
    if not T:
        raise Exception("Could not find time information")
    earliest = min(T)
    print("length <- %.2f, earliest <- %.2f" % (length, earliest), end=" ")
    R = [ math.sqrt(2*x) for x in M ]  # radius of disc
    #S = [ math.pi*2*x for x in M ]   # surface of disc
    #R = [ math.sqrt(x) for x in M ]    # radius of circle
    return T, R


def get_type(confname):
    """
    Get some value from the file
    """
    try:
        f = open(confname, 'r')
    except:
        return -2
    for line in f:
        s = line.find('preconfig.')
        if s > 0:
            s = s + len('preconfig.')
            [key, equal, val] = line[s:-1].partition('=')
            print(key, val, end=" ")
            if key == 'type':
                return int(val)
            elif key == 'mod':
                return int(val)
            else:
                print("  Warning: config.cym has " + line)
    f.close()
    return -1


def process(dirpath, momfile):
    """
        Process given directory
    """
    os.chdir(dirpath)
    if not os.path.isfile(momfile) and os.path.isfile('properties.cmp'):
        args = ['reportN', 'fiber:moment,couple,single']
        subprocess.call(args, stdout=open(momfile, 'w'))
    try:
        f = open(momfile, 'r')
    except IOError as e:
        sys.stderr.write("Error in `%s`: %s\n"%(dirpath, str(e)))
        return [dirpath, 0, 0, 0, 0, 0];
    time, data = get_radius(f)
    if do_plot:
        nice_plot(time, data)
        plt.title(dirpath, fontsize=fts)
        plt.savefig('size.png', dpi=100)
        #plt.show()
        plt.close()
    if add_fit:
        (A, B) = fit_curve(time, data)
        C = A*math.exp(-earliest*B)
    else:
        A = 0
        B = 0
        C = 0
    D = get_type('config.cym')
    return [dirpath, A, B, C, D, length];


def print_results(data, filename=''):
    """
        Save numeric data from file
    """
    f = sys.stdout
    if filename:
        f = open(filename, 'w')
    f.write("% path fit_size fit_rate size0 type polymer\n")
    for i in data:
        f.write("%12s %12.3f %12f %12f  %i %12.2f\n" % (i[0], i[1], i[2], i[3], i[4], i[5]))
    if filename:
        f.close()

#-------------------------------------------------------------------------------

def main(args):
    paths = []
    output = 'rates.txt'
    for arg in args:
        if os.path.isdir(arg):
            paths.append(arg)
        else:
            sys.stderr.write("  Warning: unexpected argument `%s'\n" % arg)
    if not paths:
        sys.stderr.write("Please specify some directory paths\n")
        sys.exit()
    cdir = os.getcwd()
    output = os.path.join(cdir,output)
    results = []
    for p in paths:
        sys.stdout.write("\n"+p+": ")
        res = process(p, 'mom.txt')
        results.append(res)
        os.chdir(cdir)
    if not os.path.isfile(output):
        print_results(results, output)


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1].endswith("help"):
        print(__doc__)
    else:
        main(sys.argv[1:])

