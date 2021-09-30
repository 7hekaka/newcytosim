#!/usr/bin/env python
#
# A script to submit analysis jobs to the SLURM queuing system
#
# Derived from submit_slurm.py
# F. Nedelec, 4.11.2020

"""
    Submit a job to the SLURM system to be called in multiple directories
    
Syntax:
    
    submit_one.py ARG [mem=????] [queue=????] [hours=INT] [days=INT] dir1 [dir2] [dir3] [...]
    
    The amount of requested memory (default=2G) should be specified in MB:
       mem=1024 (for 1 GB)
       mem=512  (for 512 MB)
       ...
    
Example:
    
    submit_one.py 'report platelet > platelet.txt' run????
    Submit one job to run command in the directories provided
    
F. Nedelec, Last updated 30.09.2021
"""


import sys, os, subprocess, tempfile

# default parameters for submission:
submit  = 'sbatch'
queue   = 'skylake'
runtime = '1:00:00'   # 1 hour
memory  = '4096'      # in MB
ncpu    = 1           # nb of threads per job

# where output is sent:
out = sys.stderr

#-------------------------------------------------------------------------------

def execute(cmd):
    """execute given command with subprocess.call()"""
    try:
        val = subprocess.call(cmd)
        if val:
            out.write("ERROR: command failed with value %i\n" % val)
            print(cmd)
    except OSError:
        out.write("ERROR, command failed: "+' '.join(cmd)+"\n")


def write_script(fd, cmd):
    """create an executable file containing the commands"""
    fid = os.fdopen(fd, "w")
    fid.write("#!/bin/bash\n")
    for s in cmd:
        fid.write(s+'\n')
    fid.close()


def sub(file):
    """return command that will submit one job"""
    # specify memory, shell, minimum number of cores and queue
    cmd  = [submit, '--nodes=1', '--ntasks=1']
    # specify number of threads if executable is threaded:
    if ncpu > 1:
        cmd += ['--cpus-per-task=%i' % ncpu]
    cmd += ['--partition='+queue]
    cmd += ['--time='+runtime] 
    cmd += ['--mem='+memory]
    # define signals sent if time is exceeded:
    cmd += ['--signal=15@120']
    cmd += ['--signal=2@60']
    # request special hardware:
    cmd += ['--constraint=avx2']
    # redirect stderr and sdtout to files:
    cmd += ['--output='+file+'.out']
    cmd += ['--error='+file+'.err']
    # call script:
    cmd += [file]
    #cmd += ['rm '+file]
    return cmd

#-------------------------------------------------------------------------------

def main(args):
    """submit jobs, depending on the arguments provided"""
    global submit, memory, runtime, queue, ncpu
    
    #find submit command:
    proc = subprocess.Popen(['which', submit], stdout=subprocess.PIPE)
    if proc.wait():
        out.write("Error: submit command `"+submit+"' not found!\n")
    else:
        submit = proc.stdout.readline().strip()

    # first argument is used for go_sim.py:
    cmd = args.pop(0)

    # create job script file:
    if not os.path.isdir('log'):
        os.mkdir('log')

    job = []
    cwd = os.getcwd()
    for arg in args:
        if os.path.isdir(arg) and os.access(arg, os.X_OK):
            job += ['cd '+os.path.abspath(arg)+' && '+cmd+';']
        else:
            [key, equal, val] = arg.partition('=')
            if key == 'mem' or key == 'memory':
                memory = val
            elif key == 'cpu' or key == 'ncpu':
                ncpu = val
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
            else:
                out.write("Error: I do not understand argument `%s'\n" % arg)
                sys.exit()
    
    if memory < 128:
        out.write("Error: requested memory (%s MB) seems too low\n" % memory)
        sys.exit()

    if ncpu < 1:
        out.write("Error: number of cpu/job must be >= 1\n")
        sys.exit()

    if job:
        fd, file = tempfile.mkstemp('', '', 'log', True)
        out.write("script %s : " % file)
        write_script(fd, job)
        os.chmod(file, 0700)
        execute(sub(file))
    else:
        out.write("Error: you need to specify at least one directory\n")


#-------------------------------------------------------------------------------

if __name__ == "__main__":
    if len(sys.argv) < 2 or sys.argv[1].endswith("help"):
        print(__doc__)
    else:
        main(sys.argv[1:])

