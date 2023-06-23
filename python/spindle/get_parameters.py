#!/usr/bin/env python3
#
# get parameter from the config file
#
# F. Nedelec, Cambridge 23.06.2023


"""
Description:
    get parameter value from the config file
    
Syntax:
    get_parameters.py DIRECTORY_PATH

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

def process(dirpath):
    """
        get parameter
    """
    res = 0
    try:
        proc = subprocess.Popen(['grep', 'preconfig', 'config.cym'], stdout=subprocess.PIPE)
        data = uncode(proc.stdout.readline()).split()[0]
        #print(data)
        proc.stdout.close()
        s = data.split('=')
        res = s[1]
    except:
        res = 0
    return res


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
            sys.stderr.write("  Error: paths must be specified\n" % arg)
            sys.exit()
    cdir = os.getcwd()
    for p in paths:
        os.chdir(p)
        res = process(p)
        print(p, res)
        os.chdir(cdir)


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1].endswith("help"):
        print(__doc__)
    else:
        main(sys.argv[1:])

