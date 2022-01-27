#!/usr/bin/env python3
# A script to run simulations sequentially.
# Copyright F. Nedelec, 16.03.2016


"""
    Run simulations sequentially and analyze the results on the fly.
 
Syntax:

    screen.py executable config_file [repeat]
    
    Bracketted arguments are optional.
    [repeat] is an integer specifying the number of run for each config file.
    
    This will use 'preconfig.py' to vary the parameters, and the template 
    file should be written accordingly (see `preconfig.py' documentation)
 
F. Nedelec, 14.03.2016
"""

# Loading modules 
import os, sys, subprocess
from pyned import find_differences, uncode
import pre_config

out = sys.stdout  #open("screen.txt", 'w')
err = sys.stderr

def uncode(arg):
    try:
        if isinstance(arg, unicode):
            return str(arg.decode('utf-8'))
    except:
        pass
    return arg

#------------------------------------------------------------------------

def job(simex, config):
    # Vary parameters and generate files in folder 'cym':
    confs = pre_config.parse(config, {}, 1, 'cym')
    for conf in confs:
        # Identify differences in the configuration file:
        try:
            dif = find_differences('config.cym', path+'/config.cym')
            for x in dif:
                out.write(' %9i' % x)
        except Exception as e:
            out.write(repr(e))
        # Run simulation
        sub = subprocess.Popen([simex, '-', conf], stdout=subprocess.PIPE)
        # Get results from standard output:
        for line in sub.stdout:
            spl = uncode(line).split()
            if len(spl) < 2:
                pass
            elif spl[0] == '%':
                if spl[1] == 'start':
                    cnt = [0, 0, 0, 0, 0, 0]
                elif spl[1] == 'end':
                    out.write(" %9i %9i %9i %9i %9i %9i" % tuple(cnt))
            elif spl[0] == 'motor':
                cnt[0] = int(spl[3])
                cnt[1] = int(spl[4]) + int(spl[5])
                cnt[2] = int(spl[6])
            elif spl[0] == 'crosslinker':
                cnt[3] = int(spl[3])
                cnt[4] = int(spl[4]) + int(spl[5])
                cnt[5] = int(spl[6])
        sub.stdout.close()
        out.write("\n")
    out.flush()


#------------------------------------------------------------------------

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
                err.write("Error: executable `%s' was already specified\n" % simex)
                sys.exit()
            simex = os.path.abspath(os.path.expanduser(arg))
        elif os.path.isfile(arg):
            if config:
                err.write("Error: config `%s' was already specified\n" % config)
                sys.exit()
            config = arg
        else:
            err.write("Error: unexpected argument `%s'\n" % arg)
            sys.exit()
    
    if not executable(simex):
        err.write("Error: executable '%s' could not be found\n" % simex)
        sys.exit()

    if not config:
        err.write("You should specify a config file on the command line\n")
        sys.exit()

    try:
        os.mkdir('cym')
    except:
        pass

    #run the simulations
    for i in range(repeat):
        job(simex, config)


#------------------------------------------------------------------------


if __name__ == "__main__":
    if len(sys.argv) < 2 or sys.argv[1].endswith("help"):
        print(__doc__)
    else:
        main(sys.argv[1:])


