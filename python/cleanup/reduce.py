#!/usr/bin/env python3
#
# reduce.py
#
# copyright F. Nedelec, 2010-2013

"""
reduce.py:
    reduce the size of cytosim file 'objects.cmo' in selected directories
    
Usage:
    reduce.py [min=INTEGER] [commit] DIRECTORIES
    
The word 'commit' must be given otherwise the directories will not be processed
Specifying 'min' sets the minimum number of frames in 'objects.cmo'

Each pass reduces the number of frames by a factor 10.
A file 'reduced' is created in each processed directory.
If such a file already exists, the directory is not processed.
This is a safety mechanism to ensure each directory is processed once only.

F. Nedelec, 28.04.2011 - 20.12.2013 - 11.12.2019
"""

import sys, os, subprocess

commit = 0
min_frames = 2

err = sys.stderr
park = 'reduced'

#------------------------------------------------------------------------

def nbFrames(file='objects.cmo'):
    """count number of frame using frametool"""
    try:
        child = subprocess.Popen(['frametool', file], stdout=subprocess.PIPE)
        str = child.stdout.readline().strip().split(' ')
        child.stdout.close()
        return int(str[1])
    except:
        return 0


def process(path):
    """cleanup directory"""
    os.chdir(path)
    files = os.listdir('.')

    print('- '*32+path)
    if 'reduced' in files:
        return
    if commit:
        if 'reduced.cmo' in files:
            stat = os.stat('reduced.cmo')
            if stat.st_size > 10000:
                os.remove('objects.cmo')
                os.rename('reduced.cmo', 'objects.cmo')
                subprocess.call(['touch', 'reduced'])
                if 'messagesR.cmo' in files:
                    os.rename('messagesR.cmo', 'messages.cmo')
                print('reduced.cmo > objects.cmo')
    else:
        if 'objects.cmo' in files:
            cnt = nbFrames('objects.cmo')
            if cnt > min_frames:
                subprocess.call(['frametool', 'objects.cmo', '0:10:', 'reduced.cmo'])
                if 'messages.cmo' in files:
                    out = open('messagesR.cmo', 'w')
                    subprocess.call(['grep', '-v', '^F[0-9]*[123456789] ', 'messages.cmo'], stdout=out)
            else:
                print(' skip (objects.cmo has %i frames)' % cnt)

#------------------------------------------------------------------------

def main(args):
    global commit, min_frames
    paths = []
    
    for arg in args:
        if arg == 'commit':
            commit = 1
        elif arg.startswith('commit='):
            commit = int(arg[7:])
        elif arg.startswith('min='):
            min_frames = int(arg[4:])
        elif os.path.isdir(arg):
            paths.append(arg)
        elif arg.endswith('*'):
            import glob
            paths.extend(glob.glob(arg))
        else:
            err.write("ignored '%s' on command line\n" % arg)

    if not paths:
        paths.append('.')

    cdir = os.getcwd()
    for path in paths:
        process(path)
        os.chdir(cdir)

if __name__ == "__main__":
    if len(sys.argv)>1 and sys.argv[1].endswith("help"):
        print(__doc__)
    else:
        main(sys.argv[1:])

