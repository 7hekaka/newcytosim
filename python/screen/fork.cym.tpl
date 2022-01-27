/**
Scanning for numerical stability in the Fork angle
FJN 22.03.2019
*/

set simul system
{
    time_step = 0.01
    viscosity = 1
    kT = 0.00042
}

set space cell
{
    shape = sphere
}

set fiber filament
{
    rigidity = 100
    segmentation = 10
    display = ( plus_end = 10,2; forces=0.1,yellow )
}

set hand binder
{
    binding = 10, 0.05
    unbinding = 0, 3
    display = ( size=16 )
}

set single anchor
{
    hand = binder
    activity = fixed
    stiffness = 100
}

set couple fork
{
    flip = 0
    hand1 = binder
    hand2 = binder
    diffusion = 0
    stiffness = 100
    activity = fork

    torque = [[stiffness]], 1.0472   %  30d = 0.5234; 60d = 1.0472; 90d = 1.5708 
}

set couple link
{
    hand1 = binder
    hand2 = binder
    diffusion = 0
    stiffness = 1
}


new cell
{
    radius = 1.5
}

new filament
{
    orientation = 0 1 0
    position = 0 0 0
    length = 2
}

new filament
{
    direction = 1 0
    position = 0 0 0
    length = 2
}

new fork
{
    attach1 = filament1, 1.0
    attach2 = filament2, 1.0
}

new link
{
    attach1 = filament1, 2.0
    attach2 = filament2, 2.0
}

new anchor
{
    position = 0 0 0
    attach = filament1, 1.0
}

run 100 system

run system; report couple:link:fork * { verbose=0; }
run system; report couple:link:fork * { verbose=0; }
run system; report couple:link:fork * { verbose=0; }
run system; report couple:link:fork * { verbose=0; }
run system; report couple:link:fork * { verbose=0; }
run system; report couple:link:fork * { verbose=0; }
run system; report couple:link:fork * { verbose=0; }
run system; report couple:link:fork * { verbose=0; }

