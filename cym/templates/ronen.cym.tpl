% A branched contractile actin network, with Ronen Zaidel-Bar
% 22 March 2022 - post production improved syntax!
[[F0=random.randint(10,2000)]]
[[A0=random.randint(10,2000)]]
[[P0=int((F0+A0)*random.uniform(0.5, 3))]]
[[M0=random.randint(10,14000)]]
[[X=random.randint(10,7000)]]
[[mod=[0, 1, 2, 4, 5, 7, 8]]]%preconfig.mod=[[mod]]
[[P=int(2*P0) if mod&1 else P0]][[P=int(0.75*P0) if mod&8 else P0]]
[[F=int(2*F0) if mod&1 else F0]]
[[A=int(1.3333*A0) if mod&2 else A0]]
[[M=int(2*M0) if mod&4 else M0]]

set simul system 
{
    time_step = 0.001
    viscosity = 0.1
    binding_grid_step = 0.1
}

set space cell
{
    shape = sphere
}

new cell
{
    radius = 6
}

set fiber filament
{
    rigidity = 0.05
    segmentation = 0.22
    confine = inside, 1
    activity = grow
    growing_speed = 1
    total_polymer = [[P]]
}

set hand binder
{
    binding = 5, 0.02
    unbinding = 1, inf
    display = ( color=gray; )
}

set hand plus_motor
{
    binding = 5, 0.02
    unbinding = 1, inf
    
    activity = move
    max_speed = 0
    stall_force = 4
    display = ( color=green; )
}

set couple xlinker
{
    hand1 = binder
    hand2 = binder
    stiffness = 250
    diffusion = 100
}

set couple motor
{
    hand1 = plus_motor
    hand2 = plus_motor
    stiffness = 250
    diffusion = 100
}

set hand arp2
{
    binding = 1, 0.01
    unbinding = 0, inf
    display = ( color=white, gray; )
}

set hand arp3
{
    unbinding = 0, inf
    activity = nucleate
    nucleate = 0.5, filament, ( length=0.020; plus_end=grow; )
    nucleation_angle = 1.22
    display = ( color=red; )
}

set single nucleator
{
    hand = arp3
    diffusion = 5
}

set couple arp23
{
    hand1 = arp2
    hand2 = arp3
    diffusion = 5
    stiffness = 250
    activity = fork
    torque = 0.5, 1.22     % 1.22 radian is 70 degrees
    trans_activated = 1
}

new [[F]] nucleator
new [[A]] arp23

run 8000 system
{
    solve = 0
}
change filament
{
    segmentation = 0.055
}

run 2000 system
{
}

change filament
{
    growing_speed = 0, 0
}
change arp3
{
    nucleate = 0
}

new [[M]] motor
new [[X]] xlinker

call equilibrate

change plus_motor
{
    max_speed = 1
}

run 2000 system
{
    nb_frames = 20
}

