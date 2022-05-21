#!/usr/bin/env python3
# Copyright F. Nedelec, March 2016


"""
    Little theory of motor binding.
 
Syntax:

    predict_network.py config_file [more_config_files]
    
Example:

     predict_network.py run*/config.cym > theory.txt
  
F. Nedelec, 16.03.2016
"""

from __future__ import print_function

import os, sys, math, subprocess
import read_config

#------------------------------------------------------------------------

def square(x):
    return x * x;

def write(arg):
    print(arg, end='')


def uncode(arg):
    try:
        if isinstance(arg, unicode):
            return str(arg.decode('utf-8'))
    except:
        pass
    return arg


def get_connectors(path, filename):
    """
        Extract number of Bridging couples from the simulation
    """
    cdir = os.getcwd()
    try:
        os.chdir(path)
    except:
        raise Exception("could not find directory `%s` from config" % path)
    if not os.path.isfile(filename):
        sub = subprocess.call(['reportN', 'couple', 'frame=1'], stdout=open(filename, 'w'))
    sM = [0, 0, 0]
    sC = [0, 0, 0]
    # Get result from file:
    f = open(filename, "r")
    for line in f:
        spl = uncode(line).split()
        if len(spl) < 2:
            pass
        elif spl[0] == '%':
            pass
        elif spl[0] == 'motor':
            sM[0] = int(spl[3])
            sM[1] = int(spl[4])+int(spl[5])
            sM[2] = int(spl[6])
        elif spl[0] == 'crosslinker':
            sC[0] = int(spl[3])
            sC[1] = int(spl[4])+int(spl[5])
            sC[2] = int(spl[6])
    f.close()
    os.chdir(cdir)
    return sM, sC

#------------------------------------------------------------------------

space_radius = 15
space_volume = 0
nb_fiber     = 0
fiber_length = 5
nb_crossings = 0
nb_motor     = 0
nb_crosslink = 0


def predict_single(binding_rate, binding_range, unbinding_rate):
    ratio_fibs = nb_fiber * fiber_length * 2 * binding_range / space_volume
    AsF = ratio_fibs * binding_rate / unbinding_rate
    #print("AsF, ", AsF)
    popF = 1.0 / ( 1 + AsF )
    popA = AsF * popF
    return [popF, popA]


def predict_couple(binding_rate, binding_range, unbinding_rate):
    total_length = nb_fiber * fiber_length
    #surface area from which it is possible to bind to a fiber
    ratio_fibs = 2 * total_length * binding_range / space_volume;
    #surface area from which it is possible to bind to a crosspoint
    ratio_cros = 4 * math.pi * nb_crossings * square(binding_range) / space_volume;
    #print("ratio_fibs %f ratio_cross %f" % (ratio_fibs, ratio_cros))
    try:
        # ratio of rates can be undefined if unbinding_rate == 0:
        bind = binding_rate / unbinding_rate;
    except:
        ratio_cros = ( 2 * math.pi + 1 ) * nb_crossings * square(binding_range) / space_volume;
        # This is to analyse non-steady state for which unbinding_rate = 0:
        time_span = 1;
        fraction_bound = 1 - math.exp( -binding_rate * time_span )
        #print("   fraction_bound = %9.6f;" % fraction_bound)
        B2 = square(fraction_bound) * ratio_cros
        B1 = fraction_bound * ratio_fibs - B2
        return [1-B1-B2, B1, B2]
    AsF = ( ratio_fibs - ratio_cros ) * bind;
    GsF = ratio_cros * bind;
    BsG = bind / 2;
    #print("AsF %f GsF %f BsG %f " % ( AsF, GsF, BsG ))
    popF = 1.0 / ( 1 + AsF + GsF + BsG * GsF )
    popA = AsF * popF;
    popG = GsF * popF;
    popB = BsG * popG;
    return [popF, popA+popG, popB]


def predict_binding(pile, name):
    """
    Check for 'single' or 'couple' with specified name
    """
    com = read_config.get_command(pile, ['set', '*', name])
    if not com:
        raise Exception("could not find class `%s` in cytosim config" % name)
    mode = com.keys[1]
    if mode == 'couple':
        hand = com.value('hand1')
    elif mode == 'single':
        hand = com.value('hand')
    cnt = read_config.get_command(pile, ['new', mode, name]).cnt
    com = read_config.get_command(pile, ['set', 'hand', hand])
    try:
        binding_rate = com.value("binding", 0)
        binding_range = com.value("binding", 1)
    except:
        binding_rate = com.value("binding_rate")
        binding_range = com.value("binding_range")
    try:
        unbinding_rate = com.value("unbinding", 0)
    except:
        unbinding_rate = com.value("unbinding_rate")
    #print(name+": ", cnt, binding_rate, binding_range, unbinding_rate)
    if mode == 'couple':
        [popF, popA, popB] = predict_couple(binding_rate, binding_range, unbinding_rate)
        return cnt, [ cnt * i for i in [popF, popA, popB] ]
    if mode == 'single':
        [popF, popA] = predict_single(binding_rate, binding_range, unbinding_rate)
        return cnt, [ cnt * i for i in [popF, popA] ]


def process(config):
    global space_radius, space_volume, nb_fiber, fiber_length, nb_crossings
    #if config.startswith('run'):
    #    write('%4s ' % config[3:7])
    pile = read_config.parse(config)
    com = read_config.get_command(pile, ['set', 'space', '*'])
    geo = com.value("geometry")
    space_radius = float(geo.split()[1])
    space_volume = math.pi * square(space_radius)
    #print("space: ",  space_radius, space_volume)
    com = read_config.get_command(pile, ['new', 'fiber', '*'])
    try:
        fiber_length = com.value("length")
        nb_fiber = com.cnt
    except:
        fiber_length = com[0].value("length")
        nb_fiber = com[0].cnt + com[1].cnt
    p0 = 1.0 / square(math.pi) - 0.0235 * fiber_length / space_radius; # we previously used 0.09;
    nb_crossings = p0 * nb_fiber * ( nb_fiber - 1 ) * square( fiber_length / space_radius )
    mesh_size = fiber_length / ( 1 + 2 * nb_crossings / nb_fiber )
    if 0:
        write("nb_fiber              = %.1f\n" % nb_fiber)
        write("fiber_length          = %.1f\n" % fiber_length)
        write("nb_crossings          = %.1f\n" % nb_crossings)
        write("nb_crossings/filament = %.1f\n" % (nb_crossings/nb_fiber))
        write("mesh_size             = %.3f\n" % mesh_size)
    nb_motor, nM = predict_binding(pile, 'motor')
    nb_crosslink, nC = predict_binding(pile, 'crosslinker')
    try:
        nb_complex, nX = predict_binding(pile, 'complex')
    except:
        nb_complex = 0;
    if 0:
        write("nb_motor     = %i  %.1f\n" % ( nb_motor, nM[2] ))
        write("nb_complex   = %i  %.1f\n" % ( nb_complex, nX[2] ))
        write("nb_crosslink = %i  %.1f\n" % ( nb_crosslink, nC[2] ))
    if 1:
        # output for 'anacouple' in which the number of fibers is varied:
        write('%10i' % nb_fiber)
        for i in nM:
            write(' %9.2f' % i)
        for i in nC:
            write(' %9.2f' % i)
    if 0:
        # extract number of bridging couples from simulation:
        [path, file] = os.path.split(config)
        sM, sC = get_connectors(path, 'couple.txt')
        for i in nM:
            write(' %9.2f' % i)
        for i in sM:
            write(' %9i' % i)
        for i in nC:
            write(' %9.2f' % i)
        for i in sC:
            write(' %9i' % i)
    if 0:
        # using the numbers from the simulation
        pM = 1.0 - math.exp( -sM[2] / nb_crossings )
        pC = 1.0 - math.exp( -sC[2] / nb_crossings )
    if 1:
        # using the predicted numbers
        pM = 1.0 - math.exp( -nM[2] / nb_crossings )
        pC = 1.0 - math.exp( -nC[2] / nb_crossings )
    if 0:
        #output for 'analyze_network'
        write(' %9i %9i' % (nb_motor, nb_crosslink))
        write(' %9.5f %9.5f' % (pM, pC))
    if 0:
        #Human friendly output:
        write(' motors       %9.1f %9.1f %9.1f   pM = %5.3f\n' % (nM[0], nM[1], nM[2], pM))
        write(' crosslinkers %9.1f %9.1f %9.1f   pC = %5.3f\n' % (nC[0], nC[1], nC[2], pC))
        #print(" probs:  motor  %9.1f crosslink %9.1f" % (pM, pC))
    if 0:
        #outdated output
        write(' %9.3f ' % (nb_crossings * pM * pC))
        write(' %9.3f ' % (nb_crossings * pM * pC * ( 1 - pC )))
        write(' %9.3f ' % (nb_crossings * pM * pC * ( 1 - pC ) * ( 1 - pC )))
    write('\n')


#------------------------------------------------------------------------

def executable(arg):
    return os.path.isfile(arg) and os.access(arg, os.X_OK)


def main(args):
    files = []
    
    # parse arguments list:
    for arg in args:
        if os.path.isfile(arg):
            files.append(arg)
        elif os.path.isdir(arg):
            files.append(arg+'/config.cym')
        else:
            sys.stderr.write("unexpected argument %s\n" % arg)
            sys.exit()
    
    if not files:
        sys.stderr.write("You must specify a config file on the command line\n")
        sys.exit()

    for f in files:
        process(f)


if __name__ == "__main__":
    if len(sys.argv) < 2 or sys.argv[1].endswith("help"):
        print(__doc__)
    else:
        main(sys.argv[1:])


