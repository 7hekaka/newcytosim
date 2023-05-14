#!/usr/bin/env python3
#
# compactify.py
#
# Copyright F. Nedelec, 2023

"""
compactify.py:
    zip and compress 'run' directories found in specified path
    
Usage:
    compactify.py DIRECTORY [DIRECTORY2] ...
    
F. Nedelec, 04.02.2023, 14.05.2023
"""

import sys, os, subprocess, shutil

err = sys.stderr

#------------------------------------------------------------------------

def process(path):
    """cleanup one directory"""
    # remove log file:
    obj = path+'/log.txt'
    if os.path.isfile(obj):
        os.remove(obj)
    # Unzip object file:
    obj = path+'/objects.cmo.gz'
    if os.path.isfile(obj):
        subprocess.call(['gunzip', obj])
    # tarzip directory:
    zip = path + '.tar.gz'
    code = subprocess.call(['tar', '-czf', zip, path])
    if code:
        err.write('tar -czf '+path+' failed!')
    else:
        shutil.move(path, os.path.join('compacted', path))
        print(path+" ----> "+zip)

#------------------------------------------------------------------------

def main(args):
    paths = []
    for arg in args:
        if os.path.isdir(arg):
            paths.append(arg)
        elif arg.endswith('*'):
            import glob
            paths.extend(glob.glob(arg))
        else:
            err.write("ignored '%s' on command line\n" % arg)
    if not paths:
        paths.append('.')
        err.write("Error: you must specify directories: scan.py COMMAND PATHS\n")
        return 2
    home = os.getcwd()
    for path in paths:
        os.mkdir('compacted')
        stuff = sorted(os.scandir(path), key=lambda e: e.name)
        for e in stuff:
            if e.is_dir():
                process(os.path.join(path, e.name))
        os.chdir(home)

#------------------------------------------------------------------------

if __name__ == "__main__":
    if len(sys.argv)>1 and sys.argv[1].endswith("help"):
        print(__doc__)
    else:
        main(sys.argv[1:])

