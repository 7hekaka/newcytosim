#!/usr/bin/env python3
#
# Plot fiber length, calculated by `report fiber:mark`
#
# F. Nedelec, Cambridge 21.06.2023


"""
Description:
    Plot fiber length, calculated by:
    report3 fiber:mark > fiber_length.txt
    
Syntax:
    plot_fiber_length.py DIRECTORY_PATH

To finish the plot:
    grep -H preconfig */config.cym > A
    sed "s/config.cym:%preconfig.augmin_source=/ /" A > augmin.txt
    cut -c 4-7 -c 9- augmin.txt > A
    cut -c 4- spindle_length.txt > L
    paste A L > AL.txt

"""

#font size:
fts = 14

import sys, os, math, subprocess
import matplotlib
import matplotlib.pyplot as plt
#matplotlib.use('SVG')


def uncode(arg):
    try:
        return arg.decode('utf-8')
    except:
        return arg


def plot_fiber_length(P, N, name):
    """
        Plot fiber length as function of frame
    """
    fig = plt.figure(figsize=(4, 3))
    X = len(P[0])
    XR = range(0, X)
    M = [0, 0, 0]
    cat = [ 'kinetochore', 'augmin', 'pole' ]
    for i in [ 2, 1, 0 ]:
        L = [ p/max(1,n) for p, n in zip(P[i], N[i]) ]
        plt.plot(XR, L, label=cat[i], linewidth=4.0)
        M[i] = max(L)
    plt.xlim(0, math.ceil(X))
    plt.ylim(0, math.ceil(max(M)))
    plt.xlabel('Frame', fontsize=fts)
    plt.ylabel(r'Length ($\mu m$)', fontsize=fts)
    plt.title('Fiber Length', fontsize=fts)
    plt.legend()
    fig.tight_layout()
    plt.savefig('fiber_length.png', dpi=75)
    #plt.show()
    plt.close()


def plot_fiber_count(N, name):
    """
        Plot fiber count as a function of frame
    """
    fig = plt.figure(figsize=(4, 3))
    X = len(N[0])
    XR = range(0, X)
    cat = [ 'kinetochore', 'augmin', 'pole' ]
    M = [0, 0, 0]
    for i in [ 2, 1, 0 ]:
        plt.plot(XR, N[i], label=cat[i], linewidth=4.0)
        M[i] = max(N[i])
    plt.xlim(0, math.ceil(X))
    plt.ylim(0, math.ceil(max(M)))
    plt.xlabel('Frame', fontsize=fts)
    plt.ylabel('Count', fontsize=fts)
    plt.title('Fiber count', fontsize=fts)
    plt.legend()
    fig.tight_layout()
    plt.savefig('fiber_count.png', dpi=75)
    #plt.show()
    plt.close()


def get_data(file):
    """
        Retreive data from file
    """
    P = [[], [], []]
    N = [[], [], []]
    for line in file:
        s = uncode(line).split()
        if len(s) < 2:
            pass
        elif s[0] == '%':
            if s[1] == "frame":
                iP = [ 0, 0, 0 ]
                iN = [ 0, 0, 0 ]
            if s[1] == "end":
                #print(iD, iN)
                for i in [ 0, 1, 2 ]:
                    P[i].append(iP[i])
                    N[i].append(iN[i])
        elif len(s) == 8:
            K = int(s[0])
            iP[K] = float(s[6])
            iN[K] = int(s[1])
    return P, N


def process(dirpath):
    """
        Process given directory
    """
    os.chdir(dirpath)
    filename = 'fiber_length.txt'
    if not os.path.isfile(filename):
        args = ['report3', 'fiber:mark']
        subprocess.call(args, stdout=open(filename, 'w'))
    res = 0
    with open(filename, 'r') as f:
        P, N = get_data(f)
        #print(D, N)
    plot_fiber_length(P, N, dirpath)
    plot_fiber_count(N, dirpath)
    # calculate mean length and count for all data:
    L = [0, 0, 0]
    C = [0, 0, 0]
    for i in [ 0, 1, 2 ]:
        L[i] = sum(P[i]) / float(sum(N[i]))
        C[i] = sum(N[i]) / len(N[i])
    print(f'{dirpath} {C[0]:.1f} {L[0]:.4f} {C[1]:.1f} {L[1]:.4f} {C[2]:.1f} {L[2]:.4f}')


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
            process(p)
            os.chdir(cdir)


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1].endswith("help"):
        print(__doc__)
    else:
        main(sys.argv[1:])

