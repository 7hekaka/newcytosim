% The parM plasmid partitionning mechanism
% F. Nedelec, January 2022

set simul system
{
    time_step = 0.001
    viscosity = 0.05
    display = ( size = 1536,512 )
    precondition = 1
}

set space cell
{
    shape = capsule
    display = ( visible=2; )
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
    growing_speed = [[round(growth,4)]]
    growing_off_speed = -0.030
    growing_force = 0.5
    shrinking_speed = -0.250
    hydrolysis_rate = [[round(hydrolysis,4)]]  % length = 2.3
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

set hand binder
{
    binding = 5, [[round(range,4)]]
    unbinding = 1, 3
    activity = rescue
    rescue_prob = [[round(rate,4)]]
}

set hand motor
{
    binding = 5, [[round(range,4)]]
    unbinding = 1, 3
    activity = move
    unloaded_speed = [[round(speed,4)]]
    stall_force = 1
}

set single target
{
    hand = nucleator
    stiffness = 100
}

set couple complex
{
    hand1 = binder
    hand2 = motor
    specificity = antiparallel
    stiffness = 2
    diffusion = 1
    fast_diffusion = 1
}

set bead plasmid
{
    confine = inside, 10
    display = ( style=3; color=0x0000FFAA; )
}

new 2 plasmid
{
    radius = 0.075
    attach = [[int(target)]] target
    position = 0.75 0 0
}

new [[int(complex)]] complex

run 50000 system
{
    nb_frames = 0
}

report plasmid:position * { verbose = 0 }

restart 2
