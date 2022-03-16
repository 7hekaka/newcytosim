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

import sys

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


def transform_file(input, output, modifs):
    try:
        i = open(input, 'r')
        o = open(output, 'w')
    except IOError as e:
        sys.stderr.write("Error: %s\n"%str(e))
        return
    cnt = 0;
    for line in i:
        cnt += 1
        if cnt in modifs:
            o.write(modifs[cnt])
        else:
            o.write(line)
    i.close()
    o.close()


def process(inx):
    i = "run%04i/config.cym" % (8*inx)
    o = "mod%04i.cym" % (8*inx)
    pam = get_parameters(i)
    print(i, pam)
    n = int(pam[126].split()[1])
    a = int(pam[127].split()[1])
    print(n, a)
    mod = {}
    mod[8] = '%preconfig.mod=32\n'
    mod[123] = '    fast_diffusion = 2\n}\n'
    mod[125] = 'new %i simplex\n' % (n+a)
    transform_file(i, o, mod)


def main(args):
    for n in range(0, 1):
        process(n)

if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1].endswith("help"):
        print(__doc__)
    else:
        main(sys.argv[1:])

