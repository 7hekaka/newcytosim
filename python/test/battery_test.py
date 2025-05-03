#!/usr/bin/env python3
#
# battery_test.py
#
# Copyright F. Nedelec, March 19 2011 --- 05.10.2024

"""
 Description:

    `battery_test.py` runs a series of .cym files
    
    This is useful to validate cytosim and discover basic errors

Live mode:

    battery_test.py bin/play live cym/*.cym

Batch mode:

    battery_test.py bin/sim cym/*.cym > test.md
    
    scan.py 'play image frame=1000' run_*

F. Nedelec, March-June 2011; 02.2013; 01.2020; 07.2021; 11.2024; 03.2025
"""

import shutil, sys, os, subprocess, time

timeout = 1200  # in seconds
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
    """run file in a new separate directory, return collected output"""
    name = os.path.split(file)[1]
    wdir = 'run_'+name.partition('.')[0];
    try:
        os.mkdir(wdir)
    except OSError:
        return f"ERROR: cannot create `{wdir}'"
    os.chdir(wdir)
    shutil.copyfile(file, 'config.cym')
    res = '# ' + name
    out = open("out.txt", 'w')
    err = open("err.txt", 'w')
    sec = time.time()
    try:
        val = subprocess.call(tool, stdout=out, stderr=err, timeout=timeout)
        if val:
            msg = f'error {val}'
        else:
            msg = 'success'
    except subprocess.TimeoutExpired as e:
        msg = str(e)
        pass
    err.close()
    out.close()
    sec = time.time() - sec
    res += f': {sec:6.2f} sec : {msg}\n'
    if verbose:
        # copy standard-error:
        res += "\t\n"
        c = 0
        with open("err.txt", 'r') as f:
            for line in f:
                res += "\t"
                res += line
                c = c + 1
        if c:
            res += "\t\n"
    return res


def worker(queue, lock):
    """
    run executable taking arguments from queue
    """
    while True:
        os.chdir(home)
        try:
            t, f = queue.get(True, 1)
            #print(' queue %s %s' % (t, f))
        except:
            break;
        res = execute(t, f, 1)
        with lock:
            print(res, flush=True)

#------------------------------------------------------------------------

def main(args):
    " run cytosim for many config files"
    err = sys.stderr
    tool = args[0].split()[0]
    if os.access(tool, os.X_OK):
        tool = os.path.abspath(tool)
    else:
        err.write("Error: you must specify an executable as first argument\n")
        sys.exit()

    njobs = 8
    files = []
    live = False
    for arg in args[1:]:
        [key, equal, val] = arg.partition('=')
        if os.path.isfile(arg):
            files.append(os.path.abspath(arg))
        elif arg=='live' or arg=='live=1':
            live = True
        elif key == 'njobs' or key == 'jobs':
            njobs = int(val)
        else:
            err.write("Ignored `"+arg+"' on the command line\n")
    
    if not files:
        print("You must specify config files!")
        sys.exit()

    njobs = min(njobs, len(files))
    
    if live:
        for f in files:
            live(tool, f)
    elif njobs > 1:
        try:
            from multiprocessing import Process, Queue, Lock
        except ImportError:
            err.write("Warning: multiprocessing module unavailable\n")
        #process in parallel with child threads:
        lock = Lock()
        queue = Queue()
        for p in files:
            queue.put((tool, p))
        jobs = []
        for n in range(njobs):
            j = Process(target=worker, args=(queue,lock,))
            jobs.append(j)
            j.start()
        # wait for completion of all jobs:
        for j in jobs:
            j.join()
            j.close()
    else:
        #process sequentially:
        for f in files:
            os.chdir(home)
            res = execute(tool, f, 1)
            print(res, flush=True)


#------------------------------------------------------------------------

if __name__ == "__main__":
    if len(sys.argv)<2 or sys.argv[1].endswith("help"):
        print(__doc__)
    else:
        main(sys.argv[1:])

