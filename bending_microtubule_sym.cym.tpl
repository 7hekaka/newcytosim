% Bending of one microtubules using clamps and a vertical force at its plus end.
% Romain PIC, L'Hay-Les-Roses 2020-03-31


[[ seg = 0.1]]% [[seg]]
[[X1=[3]]][[L=X1]]% [[L]]
[[ k1 = 0.5 ]]% [[k1]]
[[ k2 = [0.2] ]]% [[k2]]
[[ k3 = 10000 ]]% [[k3]]
[[ k4 = 5000 ]]% [[k4]]
[[force=[0.5] ]]% [[force]]
[[radius=0.02]]% [[radius]]
[[time_step=[0.000005] ]]% [[time_step]]

set simul system
{
    dim = 3
    viscosity = 1
    time_step = [[time_step]]
    precondition = 0
}

set system display
{
    style = 3;
    window_size = 1280, 720;
    point_value = 0.001;
    back_color = black;
    perspective = 2;
}

set space cell 
{
    shape = sphere
    display = ( visible = 1 )
}

new cell
{
    radius = [[L/2+2]]
}


set fiber backbone
{
    rigidity     = [[k1]]
    segmentation = [[seg]]
    binding_key  = 0
    display = ( color=green; )
    end_force= 0 [[-force]] 0, plus_end
}

set fiber filament
{
    rigidity = [[k2]]
    segmentation = [[seg]]
}

set hand binder
{
    binding = 10, 0.05
    unbinding = 0, inf
    display = ( width=5; size=12; color=red; )
}
set single clamp
{
    hand = binder
    activity = fixed
    stiffness = 10000
}

set filament display
{
    color      = blue;
    coloring   = 5;
    points     = 6, 1;
    lines      = 4, 1;
}

set tubule microtubule
{
    stiffness = [[k3]], [[k4]]
    fiber = filament
    bone = backbone
}

new 1 microtubule
{
    length = [[L]]
    position = 0 0 0
    direction = 1 0 0
}

%%% CLAMPS ON THE MICROTUBULE %%%
new clamp
{
    position= [[-L/2+1]] 0 [[radius]]
    attach=fiber4,[[1]],minus_end
}
new clamp
{
    position= [[-L/2]] 0 [[radius]]
    attach=fiber4,0,minus_end
}
new clamp
{
    position= [[-L/2+1]] 0 [[-radius]]
    attach=fiber11,[[1]],minus_end
}
new clamp
{
    position= [[-L/2]] 0 [[-radius]]
    attach=fiber11,0,minus_end
}


new clamp
{
    position= [[-L/2+1]] [[radius]] 0 
    attach=fiber14,[[1]],minus_end
}
new clamp
{
    position= [[-L/2]] [[radius]] 0 
    attach=fiber14,0,minus_end
}
new clamp
{
    position= [[-L/2+1]] [[-radius]] 0 
    attach=fiber8,[[1]],minus_end
}
new clamp
{
    position= [[-L/2]] [[-radius]] 0 
    attach=fiber8,0,minus_end
}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


run [[int(15/time_step)]] system
{
    nb_frames = 100
}

repeat 100
{
	run [[int(5/time_step/100)]] system
	{
	    nb_frames = 1
	}
	report fiber:force fiber_positions.txt
}
