#! /bin/bash
#
# This file is double-clickable on Mac-OSX, and launches cytosim/play
#
# Usage: place this file into a folder, which should contain:
#       - a cytosim display executable called 'play' or 'cytosim'
#
# Rename the file to match the folder where the simulation results are stored.
# The name 'run0001-*' will also play 'run0001' (eg. 'play0001-fast')
#
# Example: a file 'run0001.command' will launch play in 'run0001', etc.
#
#
# F. Nedelec, 2010 - April 2017


echo "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -"

############# go to directory containing the script

cd "$(dirname $0)" || exit 1;
#echo "...now in directory" $PWD

############# find executable:

exe=cytosim
candidates="play ../play play2 cytosim";
for str in $candidates; do
    if [[ -x $str ]]; then
        exe=$(pwd)/$str;
        break;
    else
        str=$(which $str);
        if [[ -x $str ]]; then
            exe=$str;
            break;
        fi
    fi
done


############# go to directory specified in the name of the script:

cmd=$(basename $0)
cmd=${cmd%.command}
cd ${cmd} || cd ${cmd%%-*} || exit 1;
#echo "...now in directory" $PWD


if [[ ! -x $exe ]]; then
    echo Error: missing executable "$exe"
    exit 1;
fi

############# find setup

args=""
candidates="display.cms ../display.cms style.cms ../style.cms";
for str in $candidates; do
    if [[ -f $str ]]; then
       args=$str;
       break;
    fi
done

############# start executable

echo Running "$exe" in "$PWD"

$exe $args;

############# clean-up:


if [ $? == 0 ]; then
    printf "\nThank you for using Cytosim!\n\n";
else
    printf "\nCytosim did not terminate normally:\n";
    printf "  Please report this problem to feedbackATcytosimDOTorg,\n"
    printf "  and attach the configuration file to the mail\n\n";
fi
