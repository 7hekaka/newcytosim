#!/usr/bin/env python3
#
# compactify.py
#
# Copyright F. Nedelec, 2023

"""
compactify.py:
    zip and compress 'run' directories specified as arguments
    
Usage:
    compactify.py run???? ...
    
F. Nedelec, 04.02.2023, 14.05.2023, 15.08.2023
"""

import sys, os, subprocess, shutil

err = sys.stderr
park = 'compacted'

#------------------------------------------------------------------------

def cleanup(path):
    """remove log files"""
    obj = [ path+'/log.txt', path+'/out.txt', path+'/err.txt', path+'/reduced' ]
    for o in obj:
        if os.path.isfile(o):
            os.remove(o)


def process(path):
    """compress one directory"""
    cleanup(path)
    # Unzip object file if present:
    obj = path+'/objects.cmo.gz'
    if os.path.isfile(obj):
        subprocess.run(['gunzip', obj], check=True)
    # tarzip directory:
    zip = path + '.tar.gz'
    code = subprocess.call(['tar', '-czf', zip, path])
    if code:
        err.write('tar -czf '+path+' failed!')
    else:
        shutil.move(path, os.path.join(park, path))
        print(path+" ----> "+zip)

#------------------------------------------------------------------------

def main(args):
    paths = []
    for arg in args:
        if os.path.isdir(arg):
            paths.append(arg)
        else:
            err.write("ignored '%s' on command line\n" % arg)
    if not paths:
        paths.append('.')
        err.write("Error: you must specify directories: compactify.py PATHS\n")
        return 2
    try:
        os.mkdir(park)
    except FileExistsError:
        pass
    paths.remove(park)
    home = os.getcwd()
    for path in paths:
        process(path)
        os.chdir(home)

#------------------------------------------------------------------------

if __name__ == "__main__":
    if len(sys.argv)>1 and sys.argv[1].endswith("help"):
        print(__doc__)
    else:
        main(sys.argv[1:])

