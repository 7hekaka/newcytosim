
# Testing Cytosim

It is good practice to be able to validate a version of the code by running a set of config files of reliable outcome.
Here we are testing the final executable.

# Current testing method

## A. Preparing the executable

Running the whole set of tests takes hours. 
 To speed up, a modification is made to 'src/sim/parser.cc', around line 700 :

	void Parser::parse_run(std::istream& is)
	{
		...
         else if ( opt.has_key("nb_steps") )
             throw InvalidSyntax("the number of steps was specified twice");
         
	+     // reduce number of steps, to speed up validation:
	+     cnt /= 16;
	+        
         if ( opt.empty() )
             execute_run(cnt);
         else
	}

This way the effective running time is reduced.

Compile the code with **assertions enabled** to make 'sim'.
This can be done with `DIM=2` or `DIM=3`.


## B. Running the test 

	mkdir test
	cd test
	cp ~/code/cytosim/bin/sim .
	cp ~/code/cytosim/python/test/battery_test.py .
	mkdir cym
	cp ~/code/cytosim/cym/*.cym cym/.
	
	./battery_test.py sim cym/*.cym &> test.md

This should run all the config files into a separate subfolder:
	
	run_actin
	run_amplify
	run_arp23
	...

Each subfolder should contain the usual cytosim output, and in addition:  

- the standard output `out.txt` 
- the standard error `err.txt`
- the log file `messages.out`

## C. Checking for reported errors

Examine `test.md` for error and warning messages.

Example of warning that can be ignored:
	
	---------------------------------------------------fiber_forces.cym
	 0.09 sec
	> Warning, a value was ignored: length = 6 (used once), 0 (used once), 0 (unused);
	>  in    23 | new cell
	>        24 | {
	>        25 |     length = 6, 0, 0
	>        26 | }


Example of an error that can be ignored:

	--------------------------------------------------endocytosis.cym
	 0.00 sec
	> Aborting since the config file specifies a different dimensionality:
	> Cytosim was compiled with DIM=2 but the config specifies `dim=3'

If the config specifies 'dim=3', and the code was compiled for 2D, it will abort at the first read for the config file.


In case of errors, the files in each subfolder `run_*` can be examined:  

- standard output `out.txt` and standard error `err.txt`
- and the log file `messages.out`


## D. Visually checking the results

2. Visually examine the resutls of the simulations. 

This can be automated using our python script `scan.py`:

For a 3D simulation, first generate an image for each subfolder:

	scan.py 'play3 image frame=10000' run*

For a 2D simulation:

	scan.py 'play2 image frame=10000' run*
 	
Possibly convert the images using [ImageMagick](https://imagemagick.org/):

	for d in run_*/image*ppm; do convert $d $d.png; done


It may help to build a HTML summary page, using our script `make_page.py`:

	make_page.py tile=6 width=256 run_*

Open and examine the results

# Todo

The script `battery_test.py` can be replaced by some automatic workflow controller.


# Idea of future development 

## Automatically evaluating the results

We could automatically compare the images, or directly the trajectory files. 

Since the simulation is stochastic, the output of every run is generally different. However, we can set the 'random_seed' to a fixed value to avoid this problem, but in addition the simulations would need to run on the same hardware, since the random number generator gives different output on little and big-endian machines.

	set simul system
	{
	    dim = 2
	    random_seed = 1 
	}

Alternatively, we can modify the source code to always initialize the `random_seed`:
In file `simul_prop.cc`:

	void SimulProp::clear()
	{
	    random_seed = 1;
	}

Possibly two series of tests need to be performed: 

- with a true random seed, as described above, to widely check for bugs
- with known seed (`random_seed = 1`), to allow for comparison with earlier results

## About this file

FJN 4.09.2020, 5.10.2024