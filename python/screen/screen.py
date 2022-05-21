#!/usr/bin/env python3
# A script to run simulations on the fly sequentially.
# Copyright F. J. Nedelec, 22.3.2019, 27.01.2022


"""
    Run simulations sequentially and analyze the results on the fly.
 
Usage:

    screen.py executable template_config_file [repeat]
    
    [repeat] is an optional integer specifying the number of run for each config file.
    
    This will use 'preconfig.py' to vary parameter values, and the template
    file should be written accordingly (see `preconfig.py' documentation)
    Copy `preconfig.py' in your working directory.

F. Nedelec, 30.01.2020, 27.01.2022
"""

try:
    import os, sys, math, random, subprocess, preconfig
except ImportError as e:
    sys.stderr.write("Error loading module: %s\n"%str(e))
    sys.exit()

out = sys.stdout  #open("screen.txt", 'w')
err = sys.stderr

def error(arg):
    err.write(arg)
    sys.exit()

#------------------------------------------------------------------------

def run_one(simex, conf):
    """
        Run simulation in current working directory and return its output
    """
    try:
        # simulation started as a subprocess
        sub = subprocess.Popen([simex, '-', conf], stdout=subprocess.PIPE)
        # get output from simulation
        out, err = sub.communicate()
        if err:
            error("Subprocess failed in %s: %s\n" %(os.getcwd(), err))
        return out.decode()
    except Exception as e:
        err.write("Failed: %s\n" % repr(e))
    return ''


def extract_data(text):
    """ This reads the output of 'report fiber:length' """
    val = []
    for line in text.split('\n'):
        s = line.split()
        if len(s) == 7 and s[0] == 'filament':
            val.append(float(s[2]))
    # return mean of the length:
    if val:
        return sum(val)/len(val)
    return 0


def run_job(simex, config, values):
    """ Generate files using preconfig and run one simulation for each file """
    files = preconfig.Preconfig().parse(config, values)
    res = []
    for conf in files:
        text = run_one(simex, conf)
        res.append(extract_data(text))
    return res


def run_many(simex, config, repeat):
    """ Vary a parameter value and calls run_job() """
    for i in range(repeat):
        X = round(random.uniform(0, 0.040), 6)
        res = run_job(simex, config, {'X':X})
        out.write("\n%9.6f    " % (X))
        for Y in res:
            out.write(" %9.4f" % Y)
        out.flush()
    out.write("\n")

#-------------------------------------------------------------------------------

def executable(arg):
    return os.path.isfile(arg) and os.access(arg, os.X_OK)

def main(args):
    config = ''
    repeat = 1
    simex = ''
    # parse arguments list:
    for arg in args:
        if arg.isdigit():
            repeat = int(arg)
        elif executable(arg):
            if simex:
                error("Error: executable `%s' was already specified\n" % simex)
            simex = os.path.abspath(os.path.expanduser(arg))
        elif os.path.isfile(arg):
            if config:
                error("Error: duplicate `%s' specified\n" % arg)
            config = arg
        else:
            error("Error: unexpected argument `%s'\n" % arg)
    # check
    if not simex:
        error("Error: executable `%s' could not be found\n" % simex)
    if not config:
        error("You must specify a template file on the command line\n")
    # run stuff
    run_many(simex, config, repeat)


#-------------------------------------------------------------------------------

if __name__ == "__main__":
    if len(sys.argv) < 2 or sys.argv[1]=='help':
        print(__doc__)
    else:
        main(sys.argv[1:])
