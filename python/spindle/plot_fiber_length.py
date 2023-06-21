#!/usr/bin/env python3
#
# Plot fiber length, calculated by `report fiber:mark`
#
# F. Nedelec, Cambridge 21.06.2023


"""
Description:
    Plot fiber length, calculated by `report fiber:mark`
    
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


def plot_data(D, N, name):
    """
        Plot surface as a function of time
    """
    fig = plt.figure(figsize=(4, 3))
    N = len(D[0])
    X = range(0, N)
    M = [0, 0, 0]
    for i in [ 0, 1, 2 ]:
        plt.plot(X, D[i], label=f"cat={i}", linewidth=4.0)
        M[i] = max(D[i])
    plt.xlim(0, math.ceil(N/100)*100)
    plt.ylim(0, math.ceil(max(M)))
    plt.xlabel('Time (s)', fontsize=fts)
    plt.ylabel('Length (um)', fontsize=fts)
    plt.title('Fiber Length', fontsize=fts)
    plt.legend()
    fig.tight_layout()


def get_data(file):
    """
        Retreive data from file
    """
    D = [[], [], []]
    N = [[], [], []]
    for line in file:
        s = uncode(line).split()
        if len(s) < 2:
            pass
        elif s[0] == '%':
            if s[1] == "frame":
                iD = [ 0, 0, 0 ]
                iN = [ 0, 0, 0 ]
            if s[1] == "end":
                #print(iD, iN)
                for i in [ 0, 1, 2 ]:
                    D[i].append(iD[i])
                    N[i].append(iN[i])
        elif len(s) == 8:
            K = int(s[0])
            iD[K] = float(s[6])
            iN[K] = int(s[1])
    return D, N


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
        D, N = get_data(f)
        #print(D, N)
        plot_data(D, N, dirpath)
        # calculate mean length for data above 1000s:
        L = [0, 0, 0]
        for i in [ 0, 1, 2 ]:
            L[i] = sum(D[i]) / float(sum(N[i]))
    plt.savefig('spindle_length.png', dpi=75)
    #plt.show()
    plt.close()
    print(f'{dirpath} {L[0]:.4f} {L[1]:.4f} {L[2]:.4f}')


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

