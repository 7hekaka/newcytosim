% A contractile actin network
% Jan 2016

set simul contract 
{
    time_step = 0.001
    viscosity = 0.1
    display = ( style=2; point_size=2; )
}

set space cell
{
    geometry = ( circle 15 )
}

new space cell

set fiber filament
{
    rigidity = 0.05
    segmentation = 0.15
    display = ( line=0.5, 2; )
}

set hand binder
{
    binding = 10, 0.01
    unbinding = 0.5, inf
    display = ( size=2; color=gray; )
}

set hand plus_motor
{
    binding = 10, 0.01
    unbinding = 0.5, inf
    
    activity = move
    max_speed = 0
    stall_force = 6
    
    display = ( size=2; color=blue; )
}

set couple crosslinker
{
    hand1 = binder
    hand2 = binder
    stiffness = 250
    diffusion = 10
    fast_diffusion = 1
}

set couple motor
{
    hand1 = plus_motor
    hand2 = plus_motor
    stiffness = 250
    diffusion = 10
    fast_diffusion = 1
}

new 2000 fiber filament
{
    length = 5
}

[[motor = random.randint(0, 10000)]]

new [[motor]] couple motor
new [[50000-5*motor]] couple crosslinker

run 2000 simul *
{
    solve = 1
}

run 18000 simul *
{
    solve = 0
}

change hand plus_motor
{
    max_speed = 0.5
}

run simul *
{   
    nb_frames = 5
    nb_steps = 5000
}
