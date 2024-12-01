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

    go_sim.py [[exe=]executable] [repeat] [preconfig=PYTHON] [park=directory] config_file [config_file]
        [config=CONFIGNAME] [argument=ARG] [argmuments=ARGS]
    
    Bracketted arguments are optional.
    If working_directory is not specified, the current directory is used.
    [repeat] is an integer specifying the number of run for each config file.
    Completed simulations will be store in the 'park' directory if specified.
    
    If a python script is specified, it should provide a function parse(input)
       You may use: preconfig.py
    
    One (or more) config files can be provided.

    Any program, even non-binary, can be run using exe=...

    The name of the config file can be specified using config=...

    Additional arguments can be passed as arg=... or args=[...]

    Examples :
        go_sim.py sim config*.cym njobs=4
        go_sim.py exe=python3.9 arg=../run_cytosim.py config*.cym

F. Nedelec, 03.2010, 10.2011, 05.2012, 04.2013, 12.2017
S. Dmitrieff 07.2022
"""

# Loading modules on the compute-cluster may fail for unexpected reasons.
# The long time.sleep() prevents zoombie nodes from accepting further jobs

try:
    import os, sys, shutil, time
except ImportError:
    host = os.getenv('HOSTNAME', 'unknown')
    sys.stderr.write("go_sim.py could not load necessary python modules on %s\n" % host)
    time.sleep(3600)
    sys.exit()

# go_sim.py ignores interupting SIGNALS, to make sure that it can properly clean-up
# if the child executable is killed on the same occasion. In this way 'CTRL-C'
# will kill the executable, but not this controlling script.

try:
    import go_sim_lib
except ImportError:
    sys.stderr.write("go_sim.py could not load go_sim_lib.py\n")
    sys.exit()

#-------------------------------------------------------------------------------

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


def is_executable(arg):
    return os.path.isfile(arg) and os.access(arg, os.X_OK)

#-------------------------------------------------------------------------------

class Gosimer:
    def __init__(self):
        #define output for error messages:
        self.err = sys.stderr
        self.out = sys.stdout

        self.repeat = 1
        self.park = ''
        self.exe = os.path.abspath('sim')
        self.script = None
        self.queue = ()
        self.arguments = []
        self.config = "config.cym"

    def run_one(self, conf, name):
        """
            run executable 'exe' with config 'conf' in a directory starting as 'name'
        """
        try:
            (val, res) = go_sim_lib.run(self.exe, conf, self.arguments, name, self.config)
            if val == 0:
                self.out.write("Completed `%s` in %s;\n" % (conf, res))
            else:
                self.out.write("Failed `%s` in %s with value %i\n" % (conf, res, val))
        except KeyboardInterrupt:
            self.err.write("go_sim.py `%s` was interrupted\n" % conf)
        # move run directory to `park` if specified:
        if os.path.isdir(self.park):
            try:
                res = go_sim_lib.park_directory(res, self.park, name)
                with open(res+"/log.txt", "a") as f:
                    f.write("parked    %s\n" % time.asctime())
            except Exception as e:
                self.err.write("go_sim.py cannot move directory: %s\n" % repr(e))

    def run_queue(self):
        """
        run items in the queue
        """
        while True:
            try:
                arg = self.queue.get(True, 1)
                #print("Queue -->", arg)
            except:
                break;
            self.run_one(*arg)

    def process(self, conf, name):
        """
            run configurations files generated from 'conf'
        """
        if not os.path.isfile(conf):
            self.err.write("go_sim.py could not find file '%s'\n" % conf)
            sys.exit()
        # generate config file(s):
        if self.script:
            import tempfile
            tmp = tempfile.mkdtemp('', '_'+name+'-', '.')
            template = tmp + '/config.cym.tpl'
            shutil.copy(conf, template)
            files = go_sim_lib.make_config(template, self.repeat, self.script, tmp)
        else:
            files = go_sim_lib.copy_config(conf, self.repeat)
        if not files:
            self.err.write("go_sim.py could not generate config files\n")
            sys.exit()
        # process all files created:
        for i, f in enumerate(files):
            if len(files) > 1:
                if name[-1].isdigit():
                    n = name + '-%04i' % i
                else:
                    n = name + '%04i' % i
            else:
                n = name
            if self.queue:
                self.queue.put((f, n))
                #print('Queued ' + f + ' ' + n)
            else:
                self.run_one(f, n)

    def main(self, args):
        root = 'run'
        files = []
        njobs = 1
        # parse arguments list:
        for arg in args:
            if arg.isdigit():
                self.repeat = int(arg)
            elif is_executable(arg):
                self.exe = os.path.abspath(arg)
            elif os.path.isfile(arg):
                files.append(arg)
            elif arg.startswith('preconfig'):
                self.script = arg
            else:
                [key, equal, val] = arg.partition('=')
                if key == 'nproc' or key == 'njobs' or key == 'jobs':
                    njobs = int(val)
                elif key == 'preconfig' or key == 'script':
                    self.script = val
                elif key == 'name':
                    root = val
                elif key == 'park':
                    self.park = val
                elif key == 'exe':
                    self.exe = val
                elif key == "args" or key == "arguments":
                    self.arguments += val
                elif key == "arg" or key == "argument":
                    self.arguments.append(val)
                elif key == "config" or key == "conf":
                    self.config = val
                else:
                    self.out.write("go_sim.py does not understand argument `%s'\n" % arg)
                    sys.exit()

        if not len(files):
            self.err.write("go_sim.py expects a config file on the command line\n")
            sys.exit()

        # prepare for multiprocessing
        if njobs > 1:
            try:
                from multiprocessing import Process, Queue
                self.queue = Queue()
            except ImportError:
                sys.stderr.write("Warning: multiprocessing unavailable\n")
        
        # process given files:
        if len(files) > 1:
            for i, f in enumerate(files):
                self.process(f, root+'%04i'%i)
        else:
            self.process(files[0], root)

        # process jobs:
        if njobs > 1:
            jobs = []
            for n in range(njobs):
                j = Process(target=self.run_queue)
                jobs.append(j)
                j.start()
            # wait for completion of all jobs:
            for j in jobs:
                j.join()
        return 0

#-------------------------------------------------------------------------------

if __name__ == "__main__":
    if len(sys.argv) < 2 or sys.argv[1].endswith("help"):
        print(__doc__)
    else:
        object = Gosimer()
        object.main(sys.argv[1:])


