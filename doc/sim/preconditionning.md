# Preconditionning

Chosing the right [precondionner](https://en.wikipedia.org/wiki/Preconditioner) can affect performance greatly.
The selection is done with the parameter 'precondition':

	set simul system
	{
	    ...
	    precondition = 0
	    verbose = 1
	}
	
Possible values of `precondition` are: 0, 1, 2, 3 and 4:

- 0 is no preconditionning
- 1 is a reduced banded symmetric preconditionner with 2 off-diagonals
- 2 is a reduced symmetric block preconditionner
- 3 is a reduced non-symmetric block preconditionner
- 4 is a fullsize non-symmetric block preconditionner

A higher setting requires more calculation/memory, while offering a better approximation of 
the true matrix, thus creating a tradeoff in terms of performance.

With `verbose=1` cytosim will report some useful information in `messages.cmo`:

	F5        947.00s   CPU      1.810s        2948s
		size 3*6293 kern 151 SMSBD +9*23660 precond 1 count  20 residual 0.00722835 
		size 3*6294 kern 151 SMSBD +9*23659 precond 1 count  16 residual 0.00813384 

where here on Frame 5: 

- `size` is the size of the linear system
- `kern` is the size of the largest mecable in the system
- `SMSBD ` is the type of matrix, in this case `Sparse Matrix Symmetric Block Diagonal`
- `precond` is the preconditionning method
- `count` is the number of step taken by the iterative solver (Biconjugate gradient stabilized)
- `residual` is the remaining error in the solution (defined by parameter `tolerance`)


# Automatic preconditionner comparison

Cytosim can automatically permutes among different preconditionner, giving a report to help the user chose the best one:

	run 3000 system
	{
	    solve = auto
	    nb_frames = 30
	}

In this case, `message.cmo` will include a summary as follows:

	F10         9.00s   CPU      0.447s           4s
	        size 3*3400 kern 17 SMSBx +9*6 precond 1 count   8 residual 0.00113671 
	        size 3*3400 kern 17 SMSBx +9*5 precond 1 count   8 residual 0.000948063
	 precond selection 8 | method count cpu | 0   34.2      4 | 1    7.8      2 | 2   12.8      5 | 4    2.0     16 |  -----> 1

Line 1 is Frame 10, at 9 seconds and the system is 3D of size `3*3400`.

- Preconditionnig 0 requires 34.2 iterations and this required ~4 sec of CPU.
- Preconditionnig 1 requires 7.8 iterations and this required ~2 sec of CPU.
- Preconditionnig 2 requires 12.8 iterations and this required ~5 sec of CPU.
- Preconditionnig 4 requires 2.0 iterations and this required ~16 sec of CPU.

Hence the best precondionning method is '1'. It allowed to solve the system about twice faster than without preconditionning (value 0).

This method is only valid if you probe the systems that is truly representative of the simulation.
One should not for instance sample only the start of the simulation, but rather its 'steady state'.


### About this file

FJN 4.09.2020
