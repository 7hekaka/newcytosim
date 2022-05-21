#!/usr/bin/env python3
#
# simplify_config.py
#
# simple extraction of data from config files
#
# Copyright F. Nedelec, 2011 - 2016

"""
    Read given config files, and print a summary of the most important values

Syntax:

    simplify_config.py FILES
    
Example:

    simplify_config run????/config.cym

Description:

    This script reads the config file, and print a subset of parameters.
    
    
F. Nedelec, Mar. 2016
"""

from __future__ import print_function

import sys, os, math
import read_config

#------------------------------------------------------------------------


def square(x):
    return x * x;


def process(filename):
    """ 
        This extracts parameters from the config file,
    """
    pile = read_config.parse(filename)
    com = read_config.get_command(pile, ['set', 'simul', '*'])
    time_step = com.value('time_step')
    com = read_config.get_command(pile, ['run', 'simul', '*'])
    try:
        com = com[-1]
    except:
        pass
    nb_frames = com.value('nb_frames')
    if not nb_frames:
        nb_frames = 1
    nb_steps = com.value('nb_steps')
    if not nb_steps:
        nb_steps = com.cnt
    #print('  time_step      = %.4f;' % time_step)
    #print('  nb_frames      = %i;' % nb_frames)
    #print('  nb_steps       = %i;' % nb_steps)
    interval = nb_steps * time_step / nb_frames;
    print('  interval       = %.2f;' % interval)
    # Here we check the important parameters for a network:
    com = read_config.get_command(pile, ['set', 'space', '*'])
    geo = com.value("geometry")
    space_radius = float(geo.split()[1])
    space_volume = math.pi * square(space_radius)
    #print('space_radius  = ', space_radius)
    com = read_config.get_command(pile, ['new', 'fiber', '*'])
    fiber_length = com.value("length")
    nb_fiber = com.cnt
    p0 = 1.0 / square(math.pi) - 0.0235 * fiber_length / space_radius; # we previously used 0.09;
    nb_crossings = p0 * nb_fiber * ( nb_fiber - 1 ) * square( fiber_length / space_radius )
    nb_crossings_per_filament = 2 * nb_crossings / nb_fiber
    min_mesh_size = fiber_length / ( 1 + nb_crossings_per_filament )
    print('  nb_fibers      = %i;' % nb_fiber, end='')
    print('  fiber_length   = %.3f;' % fiber_length)
    print('  nb_crossings   = %.2f;' % nb_crossings)
    print('  nb_crossings_per_filament = %.2f;' % nb_crossings_per_filament)
    print('  min_mesh_size  = %.5f;' % min_mesh_size)
    # Calculate buckling of filament
    com = read_config.get_command(pile, ['set', 'fiber', '*'])
    fiber_rigidity = com.value("rigidity")
    mesh_size = min_mesh_size;
    buckling_force = fiber_rigidity * 4 * square(math.pi/mesh_size)
    print('  buckling_force = %.5f;' % buckling_force, end='')
    print('  fiber_rigidity = %.5f;' % fiber_rigidity)
    com = read_config.get_command(pile, ['set', '*', 'plus_motor'])
    if com:
        stall_force = com.value("stall_force")
    else:
        stall_force = 6;
    print('  stall_force = %.2f;' % stall_force)
    beta = math.sqrt(buckling_force/stall_force)
    print('  exponent_beta = %.5f;' % beta)
 

#------------------------------------------------------------------------

def main(args):
    files = []
    
    for arg in args:
        if os.path.isdir(arg):
            files.append(arg+'/config.cym')
        elif os.path.isfile(arg):
            files.append(arg)
        else:
            sys.stderr.write("  Error: unexpected argument `%s'\n" % arg)
            sys.exit()
    
    if not files:
        files = [ 'config.cym' ]
    
    for f in files:
        if len(files) > 1:
            print('- '*32+f)
        process(f)


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1].endswith("help"):
        print(__doc__)
    else:
        main(sys.argv[1:])


