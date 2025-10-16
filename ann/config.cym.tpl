
% A classic fiber
% K. K. Acheampong, 2025
set simul system
{
    steric = 1, 500
    time_step = 0.002
    kT = 0.0042
    viscosity = 0.1 
    
}
set system display { back_color=white }

set space stripbox {
    shape = annulus
    display = ( color = blue; )
}
new stripbox {
    outer = 3.8
    inner =  2.8
    top = 8.0
    bottom = -8.0
}

%--- Crosslinker definition ---

set hand actin_binder
{
    binding_rate = 10
    binding_range = 0.060
    unbinding_rate = 0.08
    unbinding_force = 5
    display = ( color=cyan; )
}
set couple crosslinker
{
    hand1 = actin_binder
    hand2 = actin_binder
    stiffness = 10
    diffusion = 10
    specificity = parallel
     
}


%--- Filament definition ---
set fiber actin {
    rigidity = 0.075
    segmentation = 0.12
    steric          = 1, 0.025
    confine = inside, 200
    display = ( color = black; )

    
}

%--- Motor hand ---
set hand motor {
    activity = move
    binding = 10, 0.05
    unbinding = 0.1, 3
    unloaded_speed = 2.0
    stall_force = 5
    display = ( color=red; size=2 )
}

%--- Wrist motor definition ---
set single motor_on_surface {
    activity = wrist
    hand = motor
    stiffness = 10
}

[[v = [1, 10]]]

%--- Solid (vacuole) with attached motors ---
set solid vacuole1 { viscosity = [[v]]; display = ( color=blue ); steric = 1, 0.025; confine = inside, 300 }



[[num = 128]]
[[fil_len = [4, 8]]]
[[xlink_r = 10]]
[[num_fil = num if fil_len == 8 else num*2]]

new [[num_fil]] actin { length = [[fil_len]]; position = inside;  }


 
run 150000 system
{
    nb_frames = 300
}



new [[num_fil * xlink_r]] crosslinker { position = inside;  }

run 100000 system
{
    nb_frames = 200
}

[[motors = 16]]

new 20 solid vacuole1 { sphere1 = center, 0.35, [[motors]] motor_on_surface; position = inside;}



run 200000 system
{
    nb_frames = 400
}

% crosslink ratio = [[xlink_r]]
