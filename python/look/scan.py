#!/usr/bin/env python
#
# scan.py executes a given command in many directories
#
# Copyright  F. Nedelec and S. Dmitrieff 2007--2021

"""
    Execute specified command in given directories, sequentially using a given number of processes
 
Syntax:

    scan.py command directory1 [directory2] [directory3] [...] [jobs=N]
    
Example:
    
    scan.py 'play image' run* jobs=2
    
F. Nedelec, 02.2011, 09.2012, 03.2013, 01.2014, 06.2017, 07.2021
S. Dmitreff, 06.2017
"""

try:
    import sys, os, subprocess
except ImportError:
    sys.stderr.write("Error: could not load necessary python modules\n")
    sys.exit()

out = sys.stderr
verbose = 1

#------------------------------------------------------------------------

def execute(tool, path):
    """
    run executable in specified directory
    """
    os.chdir(path)
    if verbose:
        out.write('-  '*24+path+"\n")
    try:
        subprocess.call(tool, shell=True)
    except Exception as e:
        sys.stderr.write("Error: %s\n" % repr(e));


def queued(queue):
    """
    run executable taking argument from queue
    """
    while True:
        try:
            t, p = queue.get(True, 1)
        except:
            break;
        execute(t, p)


def main(args):
    """
        read command line arguments and process command
    """
    global verbose
    try:
        tool = args[0]
    except:
        out.write("Error: you should specify a command to execute\n")
        return 1

    njobs = 1
    paths = []
    for arg in args[1:]:
        if os.path.isdir(arg):
            paths.append(os.path.abspath(arg))
        elif arg.startswith('nproc=') or arg.startswith('njobs='):
            njobs = int(arg[6:])
        elif arg.startswith('jobs='):
            njobs = int(arg[5:])
        elif arg == '-':
            verbose = 0
        else:
            out.write("  Warning: unexpected argument `%s'\n" % arg)
            sys.exit()

    if not paths:
        out.write("Error: you should specify at least one directory\n")
        return 2
    
    njobs = min(njobs, len(paths))
    
    if njobs > 1:
        #process in parallel with child threads:
        try:
            from multiprocessing import Process, Queue
            queue = Queue()
            for p in paths:
                queue.put((tool, p))
            jobs = [Process(target=queued, args=(queue,)) for n in range(njobs)]
            for job in jobs:
                job.start()
            for job in jobs:
                job.join()
            # job completed:
            return 0
        except ImportError:
            out.write("Warning: multiprocessing module unavailable\n")
    #process sequentially:
    for p in paths:
        execute(tool, p)
    return 0

#------------------------------------------------------------------------

if __name__ == "__main__":
    if len(sys.argv) < 2 or sys.argv[1].endswith("help"):
        print(__doc__)
    else:
        main(sys.argv[1:])
