# Cytosim's Configuration File
 
The virtual system and the sequence of actions to be performed is specified in a configuration file (e.g. ***config.cym***). This config file must be a plain text file produced by any plain text editor.
 
Cytosim understand a [small set of commands](commands.md) and [predefined objects](objects.md) with their [associated parameters](parameters.md).
 
You will at least need to call:
 
 - `set` to define a new object category, with parameter values ([Units](units.md))
 - `new` to create objects,
 - `run` to perform simulation steps. 
 
# Example
 
Many examples can be found in the directory called ***cym***.
	 

	% Self organization of Microtubules driven by bivalent Motors
	% Adapted from Nedelec et al. Nature, 1998
	
	set simul system
	{
	    time_step = 0.01
	    viscosity = 0.05
	    display = ( style=2; )
	}
	
	set space cell
	{
	    geometry = ( circle 10 )
	}
	
	new cell
	
	set fiber microtubule
	{
	    rigidity = 20
	    segmentation = 0.5
	    confine = inside, 100
	    display = ( line_width=1; )
	}
	
	set hand kinesin
	{
	    binding_rate = 10
	    binding_range = 0.01
	    unbinding_rate = 0.1
	    unbinding_force = 3
	    
	    activity = move
	    unloaded_speed = 0.8
	    stall_force = 5
	
	    bind_also_ends = 1
	    hold_growing_end = 1
	
	    display = ( size=7; width=7 )
	}
	
	set couple complex
	{
	    hand1 = kinesin
	    hand2 = kinesin
	    stiffness = 100
	    diffusion = 10
	}
	
	new 100 microtubule
	{
	    length = 9
	}
	
	new 2000 complex
	
	set system display
	{
	    label = (Nedelec et al. 1998 -)
	}
	
	run 5000 system
	{
	    nb_frames = 50
	}
	
