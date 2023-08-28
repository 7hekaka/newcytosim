#!/usr/bin/env python3
#
# shrink.py
#
# Copyright F.J. Nedelec, 2023

"""
shrink.py:
    unzip 'run????.tar.gz' files specified, and reduce the number of frame by a factor 10
    
Usage:
    shrink.py FILE1 [FILE2] ...
    
F.J. Nedelec, 15--16.05.2023, 28.08.2023
"""

import sys, os, subprocess, shutil

err = sys.stderr
min_frames = 100

#-------------------------------------------------------------------------------

def nbFrames(file='objects.cmo'):
    """count number of frame using frametool"""
    try:
        line = subprocess.check_output(['frametool', file]).split()
        #print(line)
        return int(line[1])
    except Exception as e:
        print(e)
        return 0


def shrink(path):
    """shrink objects.cmo"""
    subprocess.run(['frametool', 'objects.cmo', '0:10:'], stdout=open('red.cmo', 'w'), check=True)
    red = nbFrames('red.cmo')
    mes = 'messages.cmo'
    if os.path.isfile(mes):
        subprocess.call(['grep', '-v', '^F[0-9]*[123456789] ', mes], stdout=open('red.mes', 'w'))
    else:
        mes = ''
    stat = os.stat('red.cmo')
    if stat.st_size > 10240 and red > 4:
        #os.rename('objects.cmo', 'full.cmo')
        os.remove('objects.cmo')
        os.rename('red.cmo', 'objects.cmo')
        if mes:
            os.rename('red.mes', mes)
        print(f': reduced {red} frames')
    else:
        print(': size(reduced.cmo) < 10kB')


def process(file, path):
    """cleanup one run"""
    cnt = 0
    try:
        # create sentinel, failing if this file already exists:
        sentinel = os.path.join(path, 'reduced')
        fd = os.open(sentinel, os.O_CREAT|os.O_EXCL|os.O_WRONLY)
    except FileExistsError as e:
        print(f'{path}: already reduced?')
        return
    print(f'{path}', end=': ')
    # Unzip object file if present:
    obj = path+'/objects.cmo.gz'
    if os.path.isfile(obj):
        subprocess.run(['gunzip', obj], check=True)
        print(f' gunzip ', end='')
    # count frames:
    cnt = nbFrames(os.path.join(path, 'objects.cmo'))
    os.write(fd, f'{cnt}\n'.encode())
    os.close(fd)
    if cnt > min_frames:
        print(f'{cnt} frames', end='')
        cdir = os.getcwd()
        os.chdir(path)
        try:
            shrink(path)
        except Exception as e:
            with open('failed', 'w') as f:
                f.write(type(e), str(e))
        os.chdir(cdir)
    else:
        print(f'{cnt} frames')


#-------------------------------------------------------------------------------

def main(args):
    global min_frames
    files = []
    paths = []
    for arg in args:
        if arg.endswith('.tarz'):
            files.append([arg, arg[:-5]])
        elif arg.endswith('.tar.gz'):
            files.append([arg, arg[:-7]])
        elif os.path.isdir(arg):
            paths.append(arg)
        elif arg.startswith('min='):
            min_frames = int(arg[4:])
        else:
            err.write("ignored '%s' on command line\n" % arg)
    if not files and not paths:
        err.write("Error: you must specify paths or files: *.tarz or *.tar.gz\n")
        return 2
    for f in files:
        try:
            # unzip file:
            subprocess.run(['tar', '-xzf', f[0]], check=True)
            print(f'{f[0]} ---> {f[1]}', end='')
            process(f[0], f[1])
        except Exception as e:
            print(type(e), e)
            err.write('tar -xzf '+file+' failed!')
            return
    for d in paths:
        process('', d)

#-------------------------------------------------------------------------------

if __name__ == "__main__":
    if len(sys.argv)>1 and sys.argv[1].endswith("help"):
        print(__doc__)
    else:
        main(sys.argv[1:])

