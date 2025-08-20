
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
    outer = 3.5
    inner =  3.0
    top = 6.0
    bottom = -6.0
}

%--- Crosslinker definition ---

set hand actin_binder
{
    binding_rate = 10
    binding_range = 0.30
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

%--- Solid (vacuole) with attached motors ---
set solid vacuole1 { viscosity = 1; display = ( color=blue ); steric = 1, 0.025; confine = inside, 300 }
set solid vacuole2 { viscosity = 2; display = ( color=green )}
set solid vacuole3 { viscosity = 0.1; display = ( color=blue ); steric = 1, 0.025; confine = inside, 300 }
set solid vacuole4 { viscosity = 0.1; display = ( color=blue ); steric = 1, 0.025; confine = inside, 300 }
set solid vacuole5 { viscosity = 0.1; display = ( color=blue ); steric = 1, 0.025; confine = inside, 300 }
set solid vacuole6 { viscosity = 0.1; display = ( color=blue ); steric = 1, 0.025; confine = inside, 300 }








new 128 actin { length = 8; position = inside;  }

 
run 50000 system
{
    nb_frames = 100
}



new 1000 crosslinker { position = inside;  }

run 50000 system
{
    nb_frames = 100
}



new 24 solid vacuole3 { sphere1 = center, 0.20, 8 motor_on_surface; position = inside;}



run 50000 system
{
    nb_frames = 100
}