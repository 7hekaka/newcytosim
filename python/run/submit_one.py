#!/usr/bin/env python3
#
# A script to submit analysis jobs to the SLURM queuing system
#
# Derived from submit.py
# F.Nedelec, 4.11.2020, 30.09.2021, 28.10.2022, 26.01.2023, 1.2.2023, 14.12.2024

"""
    Submit jobs to the SLURM system to be called in multiple directories
    One job will be submitted for each directory specified
    
Syntax:
    
    submit_one.py ARG [mem=????] [queue=????] [ncpu=INT] [hours=INT] [days=INT] dir [dir2] [...]
    
    The amount of requested memory (default=2G) should be specified in MB:
       mem=1024 (for 1 GB)
       mem=512  (for 512 MB)
       ...
    
Example:
    
    submit_one.py 'report platelet > platelet.txt' run????
    
F. Nedelec, Copyright Cambridge University. Last updated 1.2.2023
"""


import sys, os, subprocess, tempfile

# default parameters for submission:
submit  = 'sbatch'
queue   = 'icelake'
account = ''         # Project name

runtime = '3:00:00'  # 3 hours
memory  = '4096'     # in MB
ncpu    = 1          # nb of threads per job
exclusive = 0        # node exclusivity

# where error messages are sent:
err = sys.stderr

#-------------------------------------------------------------------------------

def execute(cmd):
    """execute given command with subprocess.call()"""
    try:
        val = subprocess.call(cmd)
        if val:
            err.write("ERROR: command failed with value %i\n" % val)
            print(cmd)
    except OSError:
        err.write("ERROR, command failed: "+' '.join(cmd)+"\n")


def write_script(fd, cmd):
    """create an executable file containing the commands"""
    fid = os.fdopen(fd, "w")
    fid.write("#!/bin/bash\n")
    for s in cmd:
        fid.write(s+'\n')
    fid.close()


def sub(file):
    """return script that will submit one job"""
    # specify memory, shell, minimum number of cores and queue
    cmd  = [submit, '--nodes=1', '--ntasks=1']
    # specify number of threads if executable is threaded:
    if ncpu > 1:
        cmd += ['--cpus-per-task=%i' % ncpu]
    if account:
        cmd += ['--account='+account]
    cmd += ['--partition='+queue]
    cmd += ['--time='+runtime]
    cmd += ['--mem='+memory]
    # request exclusivity on node:
    if exclusive:
        cmd += ['--exclusive']
    # define signals sent if time is exceeded:
    cmd += ['--signal=15@120']
    cmd += ['--signal=2@60']
    # redirect stderr and sdtout to files:
    cmd += ['--output='+file+'.out']
    cmd += ['--error='+file+'.err']
    # call script:
    cmd += [file]
    #cmd += ['rm '+file]
    return cmd

#-------------------------------------------------------------------------------

def executable(arg):
    return os.path.isfile(arg) and os.access(arg, os.X_OK)

def main(args):
    """submit jobs, depending on the arguments provided"""
    global submit, memory, runtime, queue, ncpu
    
    #find submit command:
    proc = subprocess.Popen(['which', submit], stdout=subprocess.PIPE)
    if proc.wait():
        err.write("Error: submit command `"+submit+"' not found!\n")
    else:
        submit = proc.stdout.readline().strip()

    # first argument is the executable:
    cmd = args.pop(0)
    if executable(cmd):
        cmd = os.path.abspath(cmd)
    
    cwd = os.getcwd()
    paths = []
    for arg in args:
        if os.path.isdir(arg) and os.access(arg, os.X_OK):
            paths.append(os.path.abspath(arg))
        else:
            [key, equal, val] = arg.partition('=')
            if key == 'mem' or key == 'memory':
                memory = val
            elif key == 'cpu' or key == 'ncpu':
                ncpu = int(val)
            elif key == 'jobs' or key == 'njobs':
                ncpu = int(val)
            elif key == 'day' or key == 'days':
                runtime = val+'-00:00:00'
            elif key == 'hour' or key == 'hours':
                runtime = val+':00:00'
            elif key == 'minute' or key == 'minutes':
                runtime = val+':00'
            elif key == 'time':
                runtime = val
            elif key == 'queue':
                queue = val
            elif key == 'account':
                account = val
            elif key == 'exclusive':
                exclusive = val
            else:
                err.write("Error: I do not understand argument `%s'\n" % arg)
                sys.exit()
    
    if int(memory) < 128:
        err.write("Error: requested memory (%s MB) seems too low\n" % memory)
        sys.exit()

    if ncpu < 1:
        err.write("Error: number of cpu/job must be >= 1\n")
        sys.exit()
    if ncpu > 128:
        err.write("Error: number of cpu/job is excessive?\n")
        sys.exit()

    for p in paths:
        fd, file = tempfile.mkstemp('', 'j', p, True)
        write_script(fd, ['cd '+p+' && '+cmd+';'])
        os.chmod(file, 0o700)
        print(file, end=' : ', flush=True)
        execute(sub(file))
    if not paths:
        err.write("Error: you need to specify at least one directory\n")


#-------------------------------------------------------------------------------

if __name__ == "__main__":
    if len(sys.argv) < 2 or sys.argv[1].endswith("help"):
        print(__doc__)
    else:
        main(sys.argv[1:])

