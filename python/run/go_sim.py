#!/usr/bin/env python3
# A script to run simulations sequentially.
# Copyright F. Nedelec, 2010--2018, Cambridge University 2019--2021
# Using multiprocessing thanks to Adolfo Alsina and Serge Dmitrieff, March 2016

"""
Synopsis:
    
    Run simulations sequentially.
    For each config file, a simulation is started in a separate 'run' directory.
    Completed runs are moved to the 'park' directory if specified.

Syntax:

    go_sim.py [executable] [repeat] [preconfig=PYTHON] [park=directory] config_file [config_file]
    
    Bracketted arguments are optional.
    If working_directory is not specified, the current directory is used.
    [repeat] is an integer specifying the number of run for each config file.
    Completed simulations will be store in the 'park' directory if specified.
    
    If a python script is specified, it should provide a function parse(input)
       You may use: preconfig.py
    
    Any number of config file can be provided (at least one).

F. Nedelec, 03.2010, 10.2011, 05.2012, 04.2013, 12.2017
"""

# Loading modules on the compute-cluster may fail for unexpected reasons.
# The long time.sleep() prevents zoombie nodes from accepting further LSF jobs

try:
    import os, sys, shutil, time
except ImportError:
    host = os.getenv('HOSTNAME', 'unknown')
    sys.stderr.write("go_sim.py could not load necessary python modules on %s\n" % host)
    time.sleep(10000)
    sys.exit()

# go_sim.py ignores interupting SIGNALS, to make sure that it can perform cleaning-up operations
# if the child executable is killed on the same occasion.
# This way 'CTRL-C' will kill the executable, but not this controlling script.

def handle_signal(sig, frame):
    sys.stderr.write("go_sim.py escaped signal %i\n" % sig)


try:
    import signal
    signal.signal(signal.SIGINT, handle_signal)
    signal.signal(signal.SIGTERM, handle_signal)
    #sys.stderr.write("go_sim.py registered its signal handler\n")
except ImportError:
    host = os.getenv('HOSTNAME', 'unknown')
    sys.stderr.write("go_sim.py could not load `signal` on %s\n" % host)
    pass


#define output for error messages:
err = sys.stderr
out = sys.stdout

njobs  = 1
repeat = 1
park   = ''
exe    = os.path.abspath('sim')

try:
    import go_sim_lib
except ImportError:
    err.write("go_sim.py could not load go_sim_lib.py\n")
    sys.exit()


#------------------------------------------------------------------------

def executable(arg):
    return os.path.isfile(arg) and os.access(arg, os.X_OK)


def run(conf, name):
    """
        run executable 'exe' with config 'conf' in a directory of name 'name'
    """
    try:
        (val, res) = go_sim_lib.run(exe, conf, name)
        if val == 0:
            out.write("Completed run `%s` in %s; " % (conf, res))
        else:
            out.write("Failed run `%s` in %s with value %i\n" % (conf, res, val))
    except KeyboardInterrupt:
        err.write("go_sim.py `%s` was interrupted\n" % conf)
    # move run directory to `park` if specified:
    if os.path.isdir(park):
        try:
            res = go_sim_lib.park_directory(res, park, name)
            with open(res+"/log.txt", "a") as f:
                f.write("parked    %s\n" % time.asctime())
        except Exception as e:
            err.write("go_sim.py cannot move directory: %s\n" % repr(e))


def run_queue(queue):
    """
    run items in the queue
    """
    while True:
        try:
            arg = queue.get(True, 1)
            run(*arg)
        except:
            break;


def process(conf, script, name, queue):
    """
        run configurations files generated from 'conf'
    """
    if not os.path.isfile(conf):
        err.write("go_sim.py could not find file '%s'\n" % conf)
        sys.exit()
    # generate config file(s):
    if script:
        import tempfile
        tmp = tempfile.mkdtemp('', '_'+name+'-', '.')
        template = tmp + '/config.cym.tpl'
        shutil.copy(conf, template)
        files = go_sim_lib.make_config(template, repeat, script, tmp)
    else:
        files = go_sim_lib.copy_config(conf, repeat)
    if not files:
        err.write("go_sim.py could not generate config files\n")
        sys.exit()
    # process all files created:
    for i, f in enumerate(files):
        if len(files) > 1:
            n = name + '-%04i' % i
        else:
            n = name
        if njobs > 1:
            queue.put((f, n))
            #print('Queued ' + f + ' ' + n)
        else:
            run(f, n)


#------------------------------------------------------------------------

def main(args):
    global njobs, repeat, park, exe, queue
    script = ''
    name   = 'run0000'
    files  = []
    njobs  = 1
    
    # parse arguments list:
    for arg in args:
        if arg.isdigit():
            repeat = int(arg)
        elif executable(arg):
            exe = os.path.abspath(arg)
        elif os.path.isfile(arg):
            files.append(arg)
        elif arg.startswith('preconfig'):
            script = arg
        else:
            [key, equal, val] = arg.partition('=')
            if key == 'nproc' or key == 'njobs' or key == 'jobs':
                njobs = val
            elif key == 'preconfig' or key == 'script':
                script = val
            elif key == 'name':
                name = val
            elif key == 'park':
                park = val
            else:
                out.write("go_sim.py does not understand argument `%s'\n" % arg)
                sys.exit()
    
    if not files:
        err.write("go_sim.py expects a config file on the command line\n")
        sys.exit()

    if not executable(exe):
        err.write("go_sim.py could not find executable '%s'\n" % exe)
        sys.exit()

    # prepare for multiprocessing
    queue = ()
    if njobs > 1:
        try:
            from multiprocessing import Process, Queue
            queue = Queue()
        except ImportError:
            out.write("Warning: multiprocessing unavailable\n")
            njobs = 1

    # process all given files:
    cnt = 0
    for conf in files:
        process(conf, script, name, queue)
        cnt += 1
        name = 'run%04i' % cnt

    # process jobs:
    if njobs > 1:
        jobs = []
        for n in range(njobs):
            j = Process(target=run_queue, args=(queue,))
            jobs.append(j)
            j.start()
        # wait for completion of all jobs:
        for j in jobs:
            j.join()
    return 0


#------------------------------------------------------------------------


if __name__ == "__main__":
    if len(sys.argv) < 2 or sys.argv[1].endswith("help"):
        print(__doc__)
    else:
        main(sys.argv[1:])


