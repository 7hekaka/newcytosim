#! /bin/bash
# This file can is double-clickable on Mac-OSX, and launches cytosim/play

####### go to directory extracted from the name of this file:

cd "${0%/*}";

####### find executable

exe=./cytosim;

if [[ ! -x $exe ]]; then
	exe=bin/play;
fi

if [[ ! -x $exe ]]; then
	echo Error: missing executable $exe
	exit 1;
fi

####### run live

$exe live;

####### salute

if [ $? == 0 ]; then
    printf "\nThank you for using Cytosim!\n\n";
else
    printf "\nCytosim did not terminate normally:\n";
    printf "  Please report this problem to feedbackATcytosimDOTorg,\n"
    printf "  and attach the configuration file to the mail\n\n";
fi
