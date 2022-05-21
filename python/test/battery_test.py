#!/usr/bin/env python3
#
# battery_test.py
#
# Copyright F. Nedelec, March 19 2011 --- 17.07.2021

"""
battery_test.py:

    run a battery of .cym files to check cytosim

Example - live:

    battery_test.py bin/play live cym/*.cym

Example - runs:

    battery_test.py bin/sim cym/*.cym
    make_image.py 'play frame=100 window_size=512,512' *_cym

F. Nedelec, March-June 2011 - Feb 2013 - Jan 2020 - July 2021
"""

import shutil, sys, os, subprocess, time

home = os.getcwd()

#------------------------------------------------------------------------

def live(tool, file):
    """run live test"""
    print(file.center(100, '~'))
    cmd = tool + ['live', file]
    val = subprocess.call(cmd)
    if val != 0:
        print('returned %i' % val)


def execute(tool, file, verbose):
    """run test in separate directory"""
    name = os.path.split(file)[1]
    wdir = 'run_'+name.partition('.')[0];
    try:
        os.mkdir(wdir)
    except OSError:
        print('skipping  '+file)
        return
    os.chdir(wdir)
    shutil.copyfile(file, 'config.cym')
    if verbose:
        print(name.rjust(100, '-'))
    out = open("out.txt", 'w')
    err = open("err.txt", 'w')
    sec = time.time()
    val = subprocess.call(tool, stdout=out, stderr=err)
    sec = time.time() - sec
    err.close()
    out.close()
    if val:
        print(' %6.2f sec : %s returned %i' % (sec, name, val))
    else:
        print(' %6.2f sec : %s' % (sec, name))
    if verbose:
        # copy standard-error:
        with open("err.txt", 'r') as f:
            for line in f:
                print("> "+line, end='')


def worker(queue):
    """
    run executable taking argument from queue
    """
    while True:
        os.chdir(home)
        try:
            t, f = queue.get(True, 1)
            #print(' queue %s %s' % (t, f))
        except:
            break;
        execute(t, f, 0)

#------------------------------------------------------------------------

def main(args):
    " run cytosim for many config files"
    tool = args[0].split()[0]
    if os.access(tool, os.X_OK):
        tool = os.path.abspath(tool)
    else:
        err.write("Error: you must specify an executable on the command line\n")
        sys.exit()

    njobs = 1
    files = []
    live = False
    err = sys.stderr
    for arg in args[1:]:
        if os.path.isfile(arg):
            files.append(os.path.abspath(arg))
        elif arg=='live' or arg=='live=1':
            live = True
        elif arg.startswith('jobs='):
            njobs = int(arg[5:])
        elif arg.startswith('njobs='):
            njobs = int(arg[6:])
        else:
            err.write("Ignored`"+arg+"' on the command line\n")
    
    if not files:
        print("You must specify config files!")
        sys.exit()

    njobs = min(njobs, len(files))
    
    if njobs > 1:
        #process in parallel with child threads:
        try:
            from multiprocessing import Process, Queue
            queue = Queue()
            for p in files:
                queue.put((tool, p))
            jobs = []
            for n in range(njobs):
                j = Process(target=worker, args=(queue,))
                jobs.append(j)
                j.start()
            # wait for completion of all jobs:
            for j in jobs:
                j.join()
            return 0
        except ImportError:
            out.write("Warning: multiprocessing module unavailable\n")
    #process sequentially:
    if live:
        for f in files:
            live(tool, f)
    else:
        for f in files:
            os.chdir(home)
            execute(tool, f, 1)


#------------------------------------------------------------------------

if __name__ == "__main__":
    if len(sys.argv)<2 or sys.argv[1].endswith("help"):
        print(__doc__)
    else:
        main(sys.argv[1:])

