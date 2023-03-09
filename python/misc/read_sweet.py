#!/usr/bin/env python3
#
# read_sweet.py
#
# read from Cytosim's packed binary files
#
# F. Nedelec, 27.08.2017, 21.11.2020


"""
Syntax:
    
    read_sweet.py FILES
    
Description:
    
    Reads Cytosim's packed binary files
    
Sweet20 format:

    [2 bytes, unsigned integer: class]
    [2 bytes, unsigned integer: info/color]
    [4 bytes, unsigned integer: identification]
    [4 bytes, float: X coordinate]
    [4 bytes, float: Y coordinate]
    [4 bytes, float: Z coordinate]

F. Nedelec, 27.08.2017, 21.11.2020
"""

import sys, os, struct


def read(path):
    """
        Read binary file with formatted data
    """
    format = 'HHIfff'
    nbytes = struct.calcsize(format) # should be 20 bytes
    #print('sweet format with %i bytes\n' % nbytes)
    with open(path, "rb") as f:
        bytes = f.read(nbytes)
        while bytes:
            data = struct.unpack(format, bytes)
            # do something here with the data:
            print(data)
            bytes = f.read(nbytes)
        f.close()


def main(args):
    paths = []
    for arg in args:
        if os.path.isfile(arg):
            paths.append(arg)
        else:
            sys.stderr.write("  Error: unexpected argument `%s'\n" % arg)
            sys.exit()
    for p in paths:
        read(p)


if __name__ == "__main__":
    if len(sys.argv) < 2 or sys.argv[1].endswith("help"):
        print(__doc__)
    else:
        main(sys.argv[1:])


