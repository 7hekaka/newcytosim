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
# F. Nedelec, 2010 - Nov. 2014


echo "- - - - - - - - - - - CYTOSIM - - - - - - - - - - - -"
echo

############# go to directory containing the script

cd $(dirname $0) || exit 1;

############# find executable:

exe=cytosim
candidates="./play ../play ../bin/play";
for str in $candidates; do
    if [[ -x $str ]]; then
        exe=$str;
        break;
    elif [[ -x $(pwd)/$str ]]; then
        exe=$(pwd)/$str;
        break;
    fi
done

if [[ ! -x $exe ]]; then
    echo Error: missing executable
    exit 1;
fi

############# go to directory specified in the name of the file:

cmd=$(basename $0)
cmd=${cmd%.command}
cd ${cmd%%-*} || cd ${cmd%%-*} || exit 1;
#echo "Now in directory" $PWD

############# find setup

args=""
candidates="display.cyp ./*.cyp ../display.cms ../*.cms";
for str in $candidates; do
    if [[ -f $str ]]; then
       args=$str;
       break;
    fi
done

############# run cytosim

$exe $args;

############# clean-up:


if [ $? == 0 ]; then
    printf "\nThank you for using Cytosim!\n\n";
else
    printf "\nCytosim did not terminate normally:\n";
    printf "  Please report this problem to feedbackATcytosimDOTorg,\n"
    printf "  and attach the configuration file to the mail\n\n";
fi
