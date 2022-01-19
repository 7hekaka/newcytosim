#!/usr/bin/env python3
#
# A script to plot for project with Ronen Zaidel-Bar
#
# F. Nedelec, Strasbourg, 13.1.2022


"""
    A script to generate additional configs for a screen

Syntax:
    
    ronen_redo.py
    
Description:
    
"""

import sys, os, math


def get_parameters(filename):
    """
        Read from config file numeric data associated with 'new'
    """
    res = {}
    cnt = 0
    with open(filename, 'r') as f:
        for line in f:
            cnt += 1
            s = line.split()
            if len(s) == 3:
                if s[0] == 'total_polymer':
                    res[cnt] = line
                elif s[0] == 'new' and s[1].isdigit():
                    res[cnt] = line
    return res


def transform_file(input, output, substitution):
    try:
        i = open(input, 'r')
        o = open(output, 'w')
    except IOError as e:
        sys.stderr.write("Error: %s\n"%str(e))
        return
    cnt = 0;
    for line in i:
        cnt += 1
        if cnt in substitution:
            o.write(substitution[cnt])
        else:
            o.write(line)
    i.close()
    o.close()
    
    
#-------------------------------------------------------------------------------


def main(args):
    for n in range(0, 128):
        fR = "run%04i/config.cym" % (7*n+1)
        ff = "run%04i/config.cym" % (7*n+5)
        d1 = "config%04i.cym" % (7*n+2)
        d2 = "config%04i.cym" % (7*n+4)
        pam = get_parameters(ff)
        pam.pop(132)
        print(ff, pam)
        mot = pam.pop(131)
        pam[8] = '%preconfig.mod=1\n'
        transform_file(fR, d1, pam)
        pam[131] = mot
        pam.pop(37)
        pam.pop(100)
        pam.pop(107)
        pam[8] = '%preconfig.mod=4\n'
        transform_file(fR, d2, pam)


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1].endswith("help"):
        print(__doc__)
    else:
        main(sys.argv[1:])

