#! /bin/bash

# This file can is double-clickable on Mac-OSX, and launches cytosim-play
#
# Usage: place this file into a folder, which should contain:
#       - a cytosim display executable called 'play' or 'cytosim'
#       - a folder 'run????', where '????' is a 4-digit number
#
# Rename the file as 'play????', matching the 'run????' folder.
#
# Example: a file 'play0001' will launch play in 'run0001', etc.
# The name 'play0001*' will also work (eg. 'play0001-nicest')
#
# F. Nedelec, May. 2010

#############

echo "- - - - - - - - - - - CYTOSIM - - - - - - - - - - - -"

############# run 'play' or 'cytosim':

exe=./cytosim;

if [[ ! -x $exe ]]; then
	exe=play;
fi

if [[ ! -x $exe ]]; then
	echo Error: missing executable $exe
	exit 1;
fi

echo " Executable is " $exe;


############# 

for cym in "$@"; do

  echo "- - - - - - - - - - - - - - - - - - " $cym
  $exe live $cym;

done



############# clean-up:

echo "Complete";
