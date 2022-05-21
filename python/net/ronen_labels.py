#!/usr/bin/env python3
#
# A script to generate lables in PNG files
#
# F. Nedelec, Cambridge, 27.2.2022, 14.03.2022


"""
    A script to generate lables in PNG files

Syntax:
    
    ronen_label.py
    
Description:
    
"""

import sys, os, math


def get_parameters(filename):
    """
        Read selected numeric data in config file:
        3668 µm polymer
        330 nucleators
        1951 arp23
        7623 motors
        1339 xlinker
    """
    res = {}
    try:
        file = open(filename, 'r')
    except IOError as e:
        sys.stderr.write("Error: %s\n"%str(e))
        return
    for line in file:
        #print(line, end='')
        s = line.split()
        if len(s) == 3:
            if s[0] == 'total_polymer':
                res['polymer'] = s[2]
            elif s[0] == 'new' and s[1].isdigit():
                res[s[2]] = s[1]
    file.close()
    return res


def make_caption_file(path):
    res = get_parameters('config.cym')
    file = open('caption.txt', 'w')
    for k, v in res.items():
        if k == 'polymer':
            file.write("%s: %s µm\\n" % (k,v))
        else:
            file.write("%s: %s\n" % (k,v))
    file.close()


def make_caption(path):
    res = get_parameters('config.cym')
    text = ''
    for k, v in res.items():
        if k == 'polymer':
            text += ". %s: %s µm\\n" % (k,v)
        else:
            text += "  %s: %s\\n" % (k,v)
    a = 'convert -size 512x512 -pointsize 48 -fill black -gravity West caption:"'
    b = '" -flatten '+path+'/caption.png'
    print(a+text+b)

#-------------------------------------------------------------------------------


def main(args):
    paths = []
    files = []
    for arg in args:
        if os.path.isdir(arg):
            paths.append(arg)
        elif os.path.isfile(arg):
            files.append(arg)
        else:
            sys.stderr.write("  Error: unexpected argument `%s'\n" % arg)
            sys.exit()
    for p in paths:
        files.append(p+'/config.cym')
    if not files:
        files = ['config.cym']
    
    if files:
        res = get_parameters(files[0])
        print("%", end=' ')
        for k, v in res.items():
            print("%10s"%k, end=' ')
        print()
        for f in files:
            res = get_parameters(f)
            for k, v in res.items():
                print("%10s"%v, end=' ')
            print()
    else:
        cwd = os.getcwd()
        for n in range(0, 50):
            path = "run%04i" % (7*n+1)
            try:
                os.chdir(path)
                make_caption(path)
                os.chdir(cwd)
            except FileNotFoundError as e:
                print(e)


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1].endswith("help"):
        print(__doc__)
    else:
        main(sys.argv[1:])

