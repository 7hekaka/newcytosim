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


def uncode(arg):
    try:
        return arg.decode('utf-8')
    except:
        return arg


def process(dirpath):
    """
        get parameter
    """
    key = []
    val = []
    filename = os.path.join(dirpath, 'config.cym')
    with open(filename, 'r') as f:
        for line in f:
            data = uncode(line)
            if data.startswith('%preconfig.'):
                s = data[11:].split('=')
                key.append(s[0])
                val.append(s[1].strip())
    return key, val


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
    keys = []
    for p in paths:
        key, val = process(p)
        if not keys:
            keys = key
            if len(paths) == 1:
                print('% path', *keys)
        if key == keys:
            print(p, *val)
        else:
            sys.stderr.write(f'  Non-uniform parameters: {key} != {keys}\n')


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1].endswith("help"):
        print(__doc__)
    else:
        main(sys.argv[1:])

