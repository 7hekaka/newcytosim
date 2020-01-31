#!/usr/bin/env python3
# A script to run simulations on the fly sequentially.
# Copyright F. J. Nedelec, 22.3.2019


"""
    Run simulations sequentially and analyze the results on the fly.
 
Syntax:

    screen.py executable template_config_file [repeat]
    
    [repeat] is an optional integer specifying the number of run for each config file.
    
    This will uses 'preconfig.py' to vary the parameters, and the template
    file should be written accordingly (see `preconfig.py' documentation)
    Copy `preconfig.py' into the directory.
 
F. Nedelec, 30.01.2020
"""

# Loading modules 
try:
    import os, sys, math, random, subprocess
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
    confs = preconfig.Preconfig().parse(config, values, repeat, '')
    # Run simulation:
    for conf in confs:
        sub = subprocess.Popen([execut, '-', conf], stdout=subprocess.PIPE)
        # Get results from standard output:
        val = []
        for data in sub.stdout:
            line = data.decode("utf-8").split()
            #print(line)
            if len(line) == 9:
                if line[1] == '1':
                    val.append(float(line[4]));
                else:
                    val.append(-float(line[4]));    
        sub.stdout.close()
        # calculate mean of the values:
        res.append(sum(val)/len(val))
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
                err.write("Error: duplicate `%s' specified\n" % arg)
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
    for k in range(1000):
        R = round(random.uniform(1,10),3)
        L = 10 #round(random.uniform(2,10),3)
        H = round(random.uniform(0,10),3)
        # rescaling factor:
        X = L * L / ( math.pi * math.pi * R );
        res = job(execut, config, {'L':L, 'H':H, 'R':R}, repeat)
        out.write("\n%9.3f %9.3f %9.3f    " % (L, H, H/L))
        for val in res:
            out.write(" %9.3f" % (val*X))
        out.flush()


#------------------------------------------------------------------------


if __name__ == "__main__":
    if len(sys.argv) < 2 or sys.argv[1]=='help':
        print(__doc__)
    else:
        main(sys.argv[1:])


