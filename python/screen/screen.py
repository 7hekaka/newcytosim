#!/usr/bin/env python3
# A script to run simulations on the fly sequentially.
# Copyright F. J. Nedelec, 22.3.2019


"""
    Run simulations sequentially and analyze the results on the fly.
 
Syntax:

    screen.py executable templat_config_file [repeat]
    
    [repeat] is an optional integer specifying the number of run for each config file.
    
    This will uses 'preconfig.py' to vary the parameters, and the template
    file should be written accordingly (see `pre_config.py' documentation)
    Copy `preconfig.py' into the directory.
 
F. Nedelec, 22.03.2019
"""

# Loading modules 
try:
    import os, sys, subprocess
except ImportError:
    sys.stderr.write("could not load essential python modules\n")
    sys.exit()

try:
    import preconfig
except ImportError:
    sys.stderr.write("could not load `preconfig'\n")
    sys.exit()


out = sys.stdout  #open("screen.txt", 'w')
err = sys.stderr

#------------------------------------------------------------------------

def job(execut, config, values, repeat):
    res = []
    # Vary parameters and generate files in current folder:
    confs = preconfig.parse(config, values, repeat)
    # Run simulation:
    for conf in confs:
        sub = subprocess.Popen([execut, '-', conf], stdout=subprocess.PIPE)
        # Get results from standard output:
        val = []
        for data in sub.stdout:
            line = data.decode("utf-8").split()
            if len(line) == 8:
                val.append(float(line[7]));
        sub.stdout.close()
        # calculate mean and variance:
        n = 0
        m = 0
        v = 0
        for x in val:
            n += 1
            m += x
            v += x * x
        m = m / n
        v = v / n - m * m
        res.append([m, v])
    return res


#------------------------------------------------------------------------

def executable(arg):
    return os.path.isfile(arg) and os.access(arg, os.X_OK)


def main(args):
    config = ''
    repeat = 1
    execut = ''
    
    # parse arguments list:
    for arg in args:
        if arg.isdigit():
            repeat = int(arg)
        elif executable(arg):
            if execut:
                err.write("Error: executable `%s' was already specified\n" % execut)
                sys.exit()
            execut = os.path.abspath(os.path.expanduser(arg))
        elif os.path.isfile(arg):
            if config:
                err.write("Error: config `%s' was already specified\n" % config)
                sys.exit()
            config = arg
        else:
            err.write("Error: unexpected argument `%s'\n" % arg)
            sys.exit()

    if not executable(execut):
        err.write("Error: executable `%s' could not be found\n" % execut)
        sys.exit()

    if not config:
        err.write("You must specify a config file on the command line\n")
        sys.exit()

    # run the simulations varying the parameter
    for k in range(100):
        dic = {'stiffness': 10*k}
        res = job(execut, config, dic, repeat)
        for val in res:
            out.write("%20s %9.3f %9.3f\n" % (k, val[0], val[1]))
        out.flush()


#------------------------------------------------------------------------


if __name__ == "__main__":
    if len(sys.argv) < 2 or sys.argv[1]=='help':
        print(__doc__)
    else:
        main(sys.argv[1:])


