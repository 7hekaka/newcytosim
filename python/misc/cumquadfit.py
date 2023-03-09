#!/usr/bin/env python3
#
# FJN, Sainsbury Laboratory, Cambridge University
# 25.9.2020


"""
    Read data and fit a quadratic model to the cumulative probability function
    
Syntax:
    
    cumquadfit.py [filename]
    
"""

import sys, math, os, random

try:
    import numpy
    import matplotlib
    import matplotlib.pyplot as plt
except ImportError:
    print("  Error: could not load matplotlib in python " + sys.version)
    sys.exit()


#set default font size for plots
font = {'family':'arial', 'weight':'normal', 'size': 18}
matplotlib.rc('font', **font)

def square(x):
    return x * x;

def random_color():
    """ Return a darkish random RGB triplet"""
    while True:
        R = random.random()
        G = random.random()
        B = random.random()
        if R + G + B < 2:
            break
    return (R, G, B)


def read(filename):
    """
        Read numeric data from file
    """
    data = []
    with open(filename, 'r') as f:
        for line in f:
            if line and line[0]!='%':
                s = list(float(x) for x in line.split())
                data.append(s)
    return data


def plot_data(X, Y):
    """ plot data on new graph """
    plt.plot(X, Y, '-', linewidth=2, color=random_color())


def plot(X, Y):
    """ plot data from multiple files and save image """
    fig = plt.figure(figsize=(9, 6))
    ax = plt.axes()
    ax.xaxis.set_major_locator(plt.MaxNLocator(12))
    ax.yaxis.set_major_locator(plt.MaxNLocator(12))
    plot_data(X, Y)
    plt.xlabel('X-axis')
    plt.ylabel('Counts')
    plt.title('Cumulative DF')
    fig.tight_layout()
    plt.show()


def sum_integer_power(N):
    sx = N * ( N + 1 ) / 2
    sx2 = N * ( N + 1 ) * ( 2*N + 1 ) / 6
    sx3 = N * N * ( N + 1 ) * ( N + 1 ) / 4
    sx4 = N * ( N + 1 ) * ( 2*N + 1 ) * ( 3*N*N + 3*N - 1 ) / 30
    return [ sx, sx2, sx3, sx4 ]


def quadratic_val(arg, pol):
    """
    return value of 2-d order polynomial
    """
    return [ pol[0] * x * x + pol[1] * x + pol[2] for x in arg ]


def quadratic_fit(data):
    """
    fit 2-d order polynomial
    """
    sx = 0; sx2 = 0; sx3 = 0; sx4 = 0;
    sx2y = 0; sxy = 0; sy = 0
    for x, y in data:
        x2 = x * x
        sx = sx + x
        sx2 = sx2 + x2
        sx3 = sx3 + x2 * x
        sx4 = sx4 + x2 * x2
        sx2y = sx2y + x2 * y
        sxy = sxy + x * y
        sy = sy + y
    mat = numpy.array([[sx4,sx3,sx2],[sx3,sx2,sx],[sx2,sx,len(data)]])
    sol = numpy.linalg.inv(mat).dot([sx2y, sxy, sy])
    return sol


def process(filename):
    """
    build cumulative distribution function and fit 2-d order polynomial
    """
    data = read(filename)
    dataX = sorted(list(zip(*data)[0]))
    dataY = [ float(y)/len(dataX) for y in range(0, len(dataX)) ]
    fit = quadratic_fit(zip(dataX, dataY))
    plot(dataX, dataY)
    plt.plot(dataX, quadratic_val(dataX, fit), '-', linewidth=1, color=(0, 0, 1))
    plt.savefig('plot', dpi=100)
    plt.close()
    print(fit)


def create(filename, cnt):
    with open(filename, 'w') as f:
        for i in xrange(cnt):
            x = random.random()
            if random.random() < x:
                f.write("%9.3f\n" % x)

#-------------------------------------------------------------------------------

def main(args):
    files = []
    name = ''
    # examine command line arguments
    for arg in args:
        if os.path.isfile(arg):
            files.append(arg)
        else:
            if not name:
                name = arg
            else:
                sys.stderr.write("Unexpected argument `%s'\n" % arg)
                sys.exit()
    # set default
    if not files:
        create(name, 8192)
    # process all files
    for f in files:
        process(f)


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1]=='help':
        print(__doc__)
    else:
        main(sys.argv[1:])

