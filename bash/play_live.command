#! /bin/bash

# This file can is double-clickable on Mac-OSX, and launches cytosim/play
#
# F. Nedelec, March. 2010

############# go the the directory containing this script:

cd ${0%/*}

echo "- - - - - - - - - - - CYTOSIM - - - - - - - - - - - -"

############# find executable:

exe=./play;

if [[ ! -x $exe ]]; then
	exe=./cytosim;
fi

if [[ ! -x $exe ]]; then
	echo Error: missing executable $exe
	exit 1;
fi

echo " Excutable is " $exe;

############# run live

$exe live;

############# salute

if [ $? == 0 ]; then
    printf "\nThank you for using Cytosim!\n\n";
else
    printf "\nCytosim did not terminate normally:\n";
    printf "  Please report this problem to feedbackATcytosimDOTorg,\n"
    printf "  and attach the configuration file to the mail\n\n";
fi
