# Executables
 
The cytosim platform has two main executables:
 
- `sim` can calculate the simulation by following instructions contained in a config file,
- `play` can generate a graphical representation of the simulation.

See also [the file types](file_types.md)

# Sim
 
`sim` performs the instructions specified in a [config file](../sim/config.md), 
and usually creates 3 output files:
 
File              |  Type  |   Content                                                  |
------------------|--------|-------------------------------------------------------------
`properties.cmp`  | text   | the attributes of the different classes of objects
`objects.cmo`     | binary | positions and states of the objects at different time points
`messages.cmo`    | text   | informations concerning the execution, such as CPU time
 
The trajectory file `objects.cmo` is by default a binary file, but it can
also be written in a text-based format (see command `write`). The property file
 `properties.cmp` is a text file containing parameters of the simulation. 
 Most of the values are copied from the config file, but parameters that specify the initial state are not included.
 
By default the config file is `config.cym` in the current working directory, but one may
give a name ending with a <b>.cym</b> on the command line to specify a different file:

	sim spindle.cym

[Check these instructions on how to use `sim`](runs.md)

# Play

`play` generates a graphical representation of the system, and it is sometimes named `cytosim`.
It can be used either to display a simulation calcualted by 'sim', or to run a live simulation.
While `sim` can run on any computer, `play` requires graphical capabilities (OpenGL),
and some of its functionalities (eg. PNG support) depend on external libraries.
 
The [display parameters](../sim/graphics.md) controlling the graphical output can
be embedded within the config file, or be specified from a second configuration file.
Most objects have a `display` parameter controlling their appearance in `play` (color, etc.).
These parameters are ignored by `sim`, but they are copied verbatim to `properties.cmp`,
from which they are read when you replay a simulation. 

The display parameters can also be specified in a separate file with extension ".cyp",
as in `play display.cyp`.
Note that the display parameters will be printed on the terminal after you press `k` in `play`,
and you can copy-paste this output to create an initial `display.cyp`.

 
`play` can be used to:

Purpose                                         |  Command                            
------------------------------------------------|-------------------------------------
Display the simulation in the current folder    |  `play`     
Display the simulation in folder `PATH`         |  `play PATH`
Start a live simulation                         |  `play live` or `play FILENAME.cym`
Generate images from a trajectory               |  `play movie`     


The live mode is automatically started if a config file extending with '.cym' is specified on the command line.

## Play - Replay mode

Use `play` to display a simulation calculated by `sim`.
By default, play reads `properties.cmp` and `objects.cmo`, in the current directory.
 
A different trajectory file can be specified:

	play final.cmo
 
## Play - Live mode
 
Use `play live` to start a simulation and display it online:

	play live
	play spindle.cym
 
This should open a window. This simulation is calculated on the fly, and nothing is saved.
Note that with `play` some [commands](../sim/commands.md) that write files to disc are disabled.

## Offscreen rendering
 
Play can generate images without opening a window, as for example:
 
Command                                |   Result                                       |
---------------------------------------|-------------------------------------------------
`play image frame=10`                  | An image representing frame 10
`play image frame=10,20,30`            | An image for each specified frames
`play image frame=10 image_format=png` | a PNG image representing frame 10 
`play image magnification=3 frame=10`  | an image for frame 10 at 3x magnified resolution
`play movie`                           | an image for each frame in the trajectory file

 
Cytosim is always able to generate <a href="http://en.wikipedia.org/wiki/Netpbm_format">PPM images</a>, 
that can be open with <a href="http://rsbweb.nih.gov/ij/">ImageJ</a>.
Cytosim may be able to generate images in the PNG formats, 
if the required library was included during compilation.
 
 
## Keyboard controls

Play is controlled by the keyboard and a pull-down menu. Here are the most useful controls:
 
Key        |   Action                                  |
-----------|--------------------------------------------
`h`        | show online help (press again to hide)
`p` `s`    | play (space bar), stop
`<` `>`    | previous, next frame
`+` `=`    | zoom in, zoom out
arrow-keys | move around
`R`        | write display parameters to terminal
`CTRL-Q`   | quit (mac osx)
 
 
# Other programs
 
 
Cytosim offers additional tools:

- [`report`](sim/report.md) extracts data from trajectory files
- [`preconfig`](https://openresearchsoftware.metajnl.com/articles/10.5334/jors.156/) generates configurations files from a template file
- `frametool` extracts frames from trajectory files

[See how they are used](runs.md).
 
# Clickable icons
 
 
On MacOS or Linux, it is possible to use a script to start cytosim with a double-click.
See for example `bash/play_live.command` and `bash/play0000.command`.

