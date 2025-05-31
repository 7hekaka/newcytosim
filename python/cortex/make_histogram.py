#!/usr/bin/env python3
# A script to scan the angle of a collision
# Maud Formanek, 2020-2021
# FJN 01.10.2021, 26.09.2022, 6.10.2022, 6.11.2022, 23.03.2023, 7.4.2023, 21.12.2024

"""
    Generate histogram of outcome of colliding microtubules,
    as a function of the angle of incidence and save as PNG image
    
    Usage: make_histogram.py [90|180] [DIRECTORY/FILE] [ANOTHER_DIRECTORY] ...
"""
    
import sys, os, random
from math import sin, cos, pi
import matplotlib.pyplot as plt
from operator import add
from swangle_lib import *
from itertools import accumulate
import statistics


do_fit = 3 # number of coefficients for the fit
do_var = 1 # add variance
do_count = 1 # add count on top of columns
n_bins = 9 # number of angular bins within [0, 90]
half_range = 1 # range is [0, 90] or [0, 180]
output = 'swangle' # name of output image file
        
# size of font used in plots:
fts = 14
# file format: `png` or `svg`
format = 'png'



def fourier_format(label, C, S):
    msg = f'Fourier {label}  COS: '
    for x in C:
        msg += f'{x:8.3f}'
    msg += '   SIN: '
    for x in S:
        msg += f' {x:8.3f}'
    return msg


def fourier_coefficients(A, Y, N):
    """ A : angle in radian, Y : periodic data with period PI """
    #A = [ a * 0.001 * pi for a in range(1000) ]
    #Y = [ sin(a) for a in A ]
    P = A[-1] + A[0] # range of data, calculated from centers of bins
    F = 2 * P / ( pi * len(Y) )
    C = [0] * N
    S = [0] * N
    for i in range(N):
        C[i] = F * sum([y*cos(2*i*a) for a, y in zip(A,Y)])
        S[i] = F * sum([y*sin(2*i*a) for a, y in zip(A,Y)])
    return C, S


def fourier_plot(ax, A, C, S, N=-1, base=[], col=(0,0,0)):
    """ ax: axis; A: angle; C, S: fourier coefficients """
    if N < 0:
        N = len(C)
    fit = [ C[0]*0.5 for a in A ]
    if base:
        fit = list(map(add, fit, base))
    for i in range(1, N):
        csi = [ C[i]*cos(2*i*a) + S[i]*sin(2*i*a) for a in A ]
        fit = list(map(add, fit, csi))
    if ax:
        D = [ a * 180 / pi for a in A ]
        ax.plot(D, fit, color=col)
    return fit


def bar_chart_deviation(counts, repeat):
    """ estimate deviation in the position of the bar chart """
    S = sum(counts)
    res = []
    if S > 0:
        V = list(range(0, len(counts)))
        P = [ x / S for x in counts ]
        for i in range(repeat):
            sel = random.choices(V, weights=P, k=S)
            N = [ sel.count(x)/S for x in V ]
            A = list(accumulate(N))
            res.append(A)
        res = zip(*res)
        res = [ statistics.stdev(x) for x in res ]
    return res;


def plot_histogram(A, rK, rX, rZ):
    """ make histogram plot of given data """
    D = [ x * 180 / pi for x in A ]
    sup = 10 * round((max(D)+D[0])/10);
    T, K, X, Z = normalize_data(A, rK, rX, rZ)
    KX = list(map(add, K, X))
    W = ( D[1] - D[0] ) * 0.94
    fig = plt.figure(figsize=(7*(2-half_range), 5))
    ax = plt.gca()
    ax.bar(D, K, W, label='K')
    ax.bar(D, X, W, bottom=K, label='X')
    ax.bar(D, Z, W, bottom=KX, label='Z')
    if 1:
        # add white horizontal bars to delineate categories:
        col = (1, 1, 1)
        L = 0.0025 # Thickness of white line separating the data
        H = [2*L]*len(A)
        ax.bar(D, H, W, bottom=[-L]*len(K), color=col)
        ax.bar(D, H, W, bottom=[x-L for x in K], color=col)
        ax.bar(D, H, W, bottom=[x-L for x in KX], color=col)
    if do_var:
        # add error bars:
        col = (1, 1, 1)
        for d, k, x, z in zip(D, rK, rX, rZ):
            t = k + x + z
            if t > 0:
                kx, xz, z = bar_chart_deviation([round(k), round(x), round(z)], 512)
                plt.errorbar(d, k/t, kx, color=col, capsize=3)
                plt.errorbar(d, (k+x)/t, xz, color=col, capsize=3)
    if do_count:
        # add count of events in each bin:
        top = ax.bar(D, [0.02]*len(A), W, bottom=1, color='white')
        for rec, cnt in zip(top, T):
            if cnt == round(cnt):
                cnt = int(cnt)
            x = rec.get_x() + rec.get_width() * 0.5
            y = rec.get_y() + rec.get_height()
            plt.text(x, y, cnt, ha='center', va='bottom')
    if do_fit:
        if half_range:
            mA = A + list(reversed([ math.pi-x for x in A ]))
            mK = K + list(reversed(K))
            mX = X + list(reversed(X))
            fineA = [ x * math.pi / 180 for x in range(0,90) ]
        else:
            mA = A
            mK = K
            mX = X
            fineA = [ x * math.pi / 180 for x in range(0,180) ]
        C, S = fourier_coefficients(mA, mK, do_fit)
        off = fourier_plot(ax, fineA, C, S, do_fit)
        print(fourier_format(' K', C, S))
        if 1:
            K90 = C[0] * 0.5
            for i in range(1, do_fit):
                K90 = K90 + C[i] * cos(i*math.pi)
            print(f' K(PI/2) = {K90:.3}')
        C, S = fourier_coefficients(mA, mX, do_fit)
        fourier_plot(ax, fineA, C, S, do_fit, off)
        print(fourier_format(' X', C, S))
    if 0:
        pink = [1,0,1]
        # 2-coefficient fit:
        C, S = fourier_coefficients(A, K, 8)
        M = C[0] - 2 * ( C[1] - C[2] + C[3] - C[4] + C[5] - C[6] + C[7] )
        N = C[2] + C[4] + C[6]
        print(M, N)
        off = fourier_plot(ax, fineA, [ M-2*N, -0.5*M, N ], [0, 0, 0], 3, [], pink)
        C, S = fourier_coefficients(A, X, 8)
        M = C[0] - 2 * ( C[1] - C[2] + C[3] - C[4] + C[5] - C[6] + C[7] )
        N = C[2] + C[4] + C[6]
        print(M, N)
        fourier_plot(ax, fineA, [ M-2*N, -0.5*M, N ], [0, 0, 0], 3, off, pink)
    if 0:
        pink = [1,0,1]
        # 2-coefficient fit with preset coefficients:
        K1 = 0.075
        K2 = 0.030
        off = fourier_plot(ax, fineA, [ 2*(K1+K2), -K1, -K2 ], [0, 0, 0], 3, [], pink)
        X1 = 0.41
        X2 = 0.12
        fourier_plot(ax, fineA, [ 2*(X1+X2), -X1, -X2 ], [0, 0, 0], 3, off, pink)
    if 0:
        # renomalize Xrossings without counting Katastrophes:
        N = [ 0 for i in A ]
        T, K, X, Z = normalize_data(A, N, rX, rZ)
        C, S = fourier_coefficients(A, X, 4)
        print(fourier_format('rX', C, S))
    # adjust axes:
    ax.set_ylim(0, 1)
    ax.set_xlim(0, sup)
    plt.xticks(fontsize=fts)
    plt.yticks(fontsize=fts)
    ax.set_ylabel('Outcome probability', fontsize=fts+2)
    ax.set_xlabel('Incoming MT angle [deg]', fontsize=fts+2)
    if format != 'svg':
        ax.legend(loc=2, fontsize=fts)
    fig.tight_layout()
    #plt.show()
    return ax


def make_histogram(data, filename, title=''):
    """ make histogram and save plot as image file """
    if half_range:
        data = fold_angles(data)
        A, K, X, Z, B = bin_collisions(data, n_bins, pi/2)
    else:
        A, K, X, Z, B = bin_collisions(data, n_bins, pi)
    events = list(zip(*data))[1]
    U = events.count('U')
    k = events.count('k')
    print_histogram(A, K, X, Z)
    print(f'({sum(B)} below, {k} early catastrophes, {U} undefined)')
    ax = plot_histogram(A, K, X, Z)
    if title:
        pass#ax.set_title(title, fontsize=fts)
    if format == 'svg':
        plt.savefig(filename+".svg", format='svg', dpi=150)
    else:
        plt.savefig(filename+".png", dpi=150)
    plt.close()


def plot_file_histogram(file):
    """plot histogram that is stored in a target file"""
    A, K, X, Z = load_histogram(file)
    print_histogram(A, K, X, Z)
    ax = plot_histogram(A, K, X, Z)
    #ax.set_title("Target", fontsize=fts)
    if format == 'svg':
        plt.savefig(file+".svg", format='svg', dpi=150)
    else:
        plt.savefig(file+".png", dpi=150)
    plt.close()
    #print("fitness = %f"%histogram_deviation(K, X, Z, target))


def main(args):
    """scan through directories"""
    global half_range, do_var, do_fit, do_count, output, format, n_bins, fts
    do_all = 0
    paths = []
    files = []
    title = ''
    for arg in args:
        [key, equal, value] = arg.partition('=')
        if arg == '180':
            half_range = 0
            n_bins = 18
            fts = 28
        elif arg == '90':
            half_range = 1
        elif arg in ('png', 'svg' ):
            format = arg
        elif key == 'histogram':
            plot_file_histogram(value)
            return
        elif key == 'title':
            title = value
        elif key == 'fit':
            do_fit = int(value)
        elif key == 'count':
            do_count = int(value)
        elif key == 'var':
            do_var = int(value)
        elif key == 'all':
            do_all = int(value)
        elif arg.endswith('.png') or arg.endswith('.svg'):
            [output, _, format] = arg.partition('.')
        elif os.path.isdir(arg):
            paths.append(arg)
        elif os.path.isfile(arg):
            files.append(arg)
        else:
            sys.stderr.write("Error: unexpected argument `%s`\n"%arg)
            return
    all = []
    cwd = os.getcwd()
    for file in files:
        if file.endswith('.hist'):
            plot_file_histogram(file)
        elif file.endswith('.txt'):
            data = read_collision(file)
            if do_all:
                make_histogram(data, file, file)
            all.extend(data)
        elif file.endswith('.csv'):
            data = read_collision_csv(file)
            if do_all:
                make_histogram(data, file, file)
            all.extend(data)
        else:
            sys.stderr.write("skipped argument `%s`\n"%file)
    for path in paths:
        os.chdir(path)
        f = 'results.txt'
        g = 'all_results.txt'
        if os.path.isfile(f):
            all.extend(read_collision(f))
        elif os.path.isfile(g):
            # will make a separate plot for this data set:
            data = read_collision(g)
            make_histogram(data, output, path)
        else:
            sys.stderr.write(f"skipped `{path}` which has no `{f}`\n")
        os.chdir(cwd)
    #print("--------------")
    if all and len(paths)+len(files):
        make_histogram(all, output, title)


#-----------------------------------------------------------------------------

if __name__ == "__main__":
    if len(sys.argv) < 2 or sys.argv[1]=='help':
        print(__doc__)
    else:
        main(sys.argv[1:])

