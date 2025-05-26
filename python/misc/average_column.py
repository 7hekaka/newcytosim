#!/usr/bin/env python3

# FJN 26.05.2025

"""
    Compute average values for columns of data

Syntax:

    average_column.py FILE [skip=INTEGER]

The option allow to skip a number of lines from the start:

        average_column.py order.txt skip=100

F. Nedelec, 26.05.2025
"""

import sys, os

def get_data(filename, skip):
    """ 
        Extracts column data from a file, skipping 'skip' lines
    """
    cnt = 0
    res = []
    with open(filename, 'r') as f:
        for L in f:
            L = L.strip('\n')
            # replace double spaces:
            while True:
                M = L.replace('  ',' ')
                if M == L:
                    break
                L = M
            # skip commented and empty lines:
            if L and L[0] != '%':
                cnt += 1
                if cnt > skip:
                    R = []
                    # convert column values to floats:
                    for S in L.split(" "):
                        try:
                            X = float(S)
                            R.append(X)
                        except:
                            R.append(S)
                    res.append(R)
        return res


def average_data(data):
    """ 
        Average data per column
    """
    res = data[0]
    cnt = len(data)
    for i, x in enumerate(res):
        if isinstance(x, float):
            #print(i,x)
            for line in data[1:]:
                res[i] += line[i]
            res[i] /= cnt
    return res
    

#-------------------------------------------------------------------------------

def main(args):
    skip = 1
    paths = []
    for arg in args:
        [key, equal, value] = arg.partition('=')
        if os.path.isfile(arg):
            paths.append(arg)
        elif key == 'skip':
            skip = int(value)
        else:
            sys.stderr.write("  Error: unexpected argument `%s'\n" % arg)
            sys.exit()
    for P in paths:
        dat = get_data(P, skip)
        res = average_data(dat)
        print(res)

 
#-------------------------------------------------------------------------------

if __name__ == "__main__":
    if len(sys.argv) < 2 or sys.argv[1].endswith("help"):
        print(__doc__)
    else:
        main(sys.argv[1:])

