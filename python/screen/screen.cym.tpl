% The parM plasmid partitionning mechanism
% F. Nedelec, January 2022

set simul system
{
    time_step = 0.001
    viscosity = 0.05
    precondition = 1
}

set space cell
{
    shape = capsule
}

new cell
{
    length = 2
    radius = 0.25
}

set fiber filament
{
    rigidity = 0.1
    segmentation = 0.1
    confine = inside, 10
    max_length = 3
    
    activity = dynamic
    growing_speed = 0.100
    growing_off_speed = [[-X]]
    growing_force = 0.5
    shrinking_speed = -0.250
    hydrolysis_rate = 0
    unhydrolyzed_prob = 0
    display = ( line_width=3; )
}

set hand nucleator
{
    binding = 0, 0.015
    unbinding = 0, 3
    activity = nucleate
    nucleate = 5, filament, ( length=0.010; end_state=grow, grow )
}

set single target
{
    hand = nucleator
    stiffness = 100
}

set bead plasmid
{
    confine = inside, 10
    display = ( style=3; color=0x0000FFAA; )
}

new 2 plasmid
{
    radius = 0.075
    attach = 3 target
    position = 0.75 0 0
}

run 50000 system
{
    nb_frames = 0
}

report fiber:length * { verbose = 0 }

