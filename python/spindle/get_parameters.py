#!/usr/bin/env python3
#
# get parameters from the config files in specified directories
#
# F. Nedelec, Cambridge 23.06.2023, Strasbourg 31.07.2023


"""
Description:
    get parameter value from the config file.
    Extract lines starting with '%preconfig.' and print values in columnar format
    
Syntax:
    get_parameters.py DIRECTORY_PATH

"""

import sys, os, subprocess


def read_metadata(filename, pattern = 'config.'):
    """
    extract lines with '%config.PARAMETER=VALUE' and
    return dictionary linking parameters to values
    """
    res = dict()
    with open(filename, 'r') as f:
        for line in f:
            s = line.find(pattern)
            if s > 0:
                s += len(pattern)
                code = line[s:-1].rstrip(';')
                [k, _, v] = code.partition('=')
                try:
                    v = float(v)
                    if v == int(v):
                        v = int(v)
                except:
                    pass
                if k:
                    res[k] = v
    return res

#------------------------------------------------------------------------

def main(args):
    paths = []
    for arg in args:
        if os.path.isdir(arg):
            paths.append(arg+'/config.cym')
        elif os.path.isfile(arg):
            paths.append(arg)
        else:
            sys.stderr.write("  Error: unexpected argument `%s'\n" % arg)
            sys.exit()
    if not paths:
            sys.stderr.write("  Error: paths must be specified\n" % arg)
            sys.exit()
    for p in paths:
        dic = read_metadata(p)
        print(p, dic)


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1].endswith("help"):
        print(__doc__)
    else:
        main(sys.argv[1:])

