# Getting numbers out of Cytosim using `report`

Cytosim has a tool called `report` to extract various information, such as coordinates of the filaments, their state or length, etc.

This program is normally compiled by `make`, and should be in the `bin` folder.
If this is not the case, try:
	
	make report

It is important to use the same settings as for `sim` and `play`, in particular to not change the DIM defined in `dim.h`.

Invoke `report` from the command line, and the output will be sent to the terminal. 
To get the coordinates of the vertices that make the fibers, invoke this command in a directory containing the output files of 'sim':

	report fiber:points
	
The output will be something like:

	%   fiber f1:100
	100     -5.88    -6.137
	100    -6.023    -5.658
	100    -6.165    -5.178
	100    -6.308    -4.699
	100    -6.453    -4.221
	100      -6.6    -3.743
	100    -6.749    -3.266
	100    -6.902     -2.79
	100     -7.06    -2.315
	100     -7.22    -1.842
	100    -7.382    -1.369
	100    -7.546   -0.8961
	100    -7.711    -0.424
	100    -7.876   0.04795
	100    -8.041    0.5199
	100    -8.207    0.9915
	100    -8.375     1.463
	100    -8.543     1.933
	100    -8.711     2.404
	% end

It is possible to redirect the output to a file:

	report fiber:points > filaments.txt

You can also restrict reporting to a particular frame:

	report fiber:points frame=10 > filaments.txt

or some selected frames:

	report fiber:points frame=10,20,30 > filaments.txt

# Framewise `report`

'reportF' is similar to 'report', but sends the report of each frame into a different file. Hence,

	reportF fiber:points

will create multiple files:

	report0000.txt
	report0001.txt
	report0002.txt
	report0003.txt
	report0004.txt
	...

There are many options to `report` and you can find a list in `src/sim/simul_report.cc`
