#!/usr/bin/env python3
# A script to run simulations sequentially.
# Copyright F. Nedelec, Jan. 2016


"""
    Run simulations sequentially and analyze the results on the fly.
 
Syntax:

    screen.py executable config_file [repeat]
    
    Bracketted arguments are optional.
    [repeat] is an integer specifying the number of run for each config file.
    
    The template file should specify parameters place-holders with a '$' sign,
    and `screen.py` needs to be updated to define the variations of the parameters
 
F. Nedelec, 01.2016
"""

# Loading modules 
import os, sys, subprocess, random, string

out = sys.stdout  #open("screen.txt", 'w')
err = sys.stderr

#------------------------------------------------------------------------

def parameters(index):
    """
        Return dictionary with varied parameter values
    """
    dic = {'off_rate'  : 1, #round(1.0 / 2**random.uniform(0, 6), 3),
           'nb_motors' : random.randint(0, 100) }
    return dic


def job(simex, template, index):
    """
        Run one simulation and ouput result as one line in 'out'
    """
    dic = parameters(index)
    # create config file:
    content = template.substitute(dic)
    with open('config.cym', 'w') as file:
        file.write(content)
        file.close()
    # run simulation:
    sub = subprocess.Popen([simex, '-', 'config.cym'], stdout=subprocess.PIPE)
    # get result from standard output:
    line = sub.stdout.readline().encode("utf-8")
    while line.startswith('%'):
        line = sub.stdout.readline().encode("utf-8")
    sub.stdout.close()
    # extract desired columns from the output line
    split = line.split()
    result = split[4] # + ' ' + split[5]
    # Print parameters and output on single line:
    for k in sorted(dic):
        out.write(str(dic[k])+"  ")
    out.write("  %s\n" % result)
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
        input = open(config, 'r').read()
        template = string.Template(input)
    except:
        err.write("could not load config file `%s'\n" % config)
        sys.exit()

    #run the simulations
    for i in range(repeat):
        job(simex, template, i)


#------------------------------------------------------------------------


if __name__ == "__main__":
    if len(sys.argv) < 2 or sys.argv[1].endswith("help"):
        print(__doc__)
    else:
        main(sys.argv[1:])


