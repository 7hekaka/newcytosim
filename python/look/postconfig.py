#!/usr/bin/env python3
#
# POSTCONFIG, extract parameter from PRECONFIG output files
#
# Copyright Francois J. Nedelec, Cambridge University 2020
# Created 28.09.2020


"""
# SYNOPSIS

    Postconfig read files created by preconfig, to extract parameters

# DESCRIPTION

    This relies on the template file printing some parameters as:

        [[ x = random.uniform(0,1) ]]
        %preconfig.x = [[ x ]];
    
    postconfig will extract the second line (eg. `%preconfig.x = 0.234`)
    and create a Python dictionnary containing {'x': 0.234}
    
"""

import os, sys

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


def main(args):
    """
        process arguments and perform corresponding task
    """
    verbose = 1
    inputs = []
    path = ''
    
    for arg in args:
        #print("postconfig argument `%s'" % arg)
        if os.path.isdir(arg):
            path = arg
        elif os.path.isfile(arg):
            inputs.append(arg)
        elif arg == '-':
            verbose = 0
        elif arg == '+':
            verbose = 2
        else:
            sys.stderr.write("  Error: unexpected argument `%s'\n" % arg)
            sys.exit()

    if not inputs:
        sys.stderr.write("  Error: you must specify an input file\n")
        sys.exit()

    for i in inputs:
        #out.write("Reading %s\n" % i)
        res = read_metadata(i)
        res['file'] = i
        print(res)
    return res


#-------------------------------------------------------------------------------

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("You must specify a file (for instructions, invoke with option '--help')")
    elif sys.argv[1].endswith("help"):
        print(__doc__)
    else:
        main(sys.argv[1:])

