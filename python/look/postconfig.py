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

import os, sys, re

PATTERN = '%preconfig.'


def extract(file, path):
    """
        Extract lines starting with PATTERN, removing '%'
    """
    res = ''
    with open(file, 'r') as f:
        for line in f:
            if line.startswith(PATTERN):
                res += line[len(PATTERN):];
    return res


def process(arg):
    """
        Read `arg` for parameter values and return dictionary
    """
    pam = {}
    #print("postconfig::execute("+arg+")")
    for str in arg.splitlines():
        try:
            res = re.match(r" *([a-zA-Z]\w*) *= *(.*)", str)
            if res and len(res.groups()) == 2:
                k = res.group(1)
                v = res.group(2).split()[0]
                try:
                    pam[k] = float(v)
                except:
                    pam[k] = v
        except Exception as e:
            sys.stderr.write("\033[95m")
            sys.stderr.write("Error evaluating `%s':\n" % arg)
            sys.stderr.write("\033[0m")
            sys.stderr.write("    "+str(e)+'\n')
            sys.exit(1)
    return pam


def main(args):
    """
        process arguments and perform corresponding task
    """
    verbose = 1
    inputs = []
    path = ''
    
    for arg in args:
        #print("preconfig argument `%s'" % arg)
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
        str = extract(i, path)
        res = process(str)
        res['file'] = i
    return res


#-------------------------------------------------------------------------------

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("You must specify a file (for instructions, invoke with option '--help')")
    elif sys.argv[1].endswith("help"):
        print(__doc__)
    else:
        res = main(sys.argv[1:])
        print(res)

