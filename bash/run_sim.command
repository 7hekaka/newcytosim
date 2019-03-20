#! /bin/bash

####### go to directory extracted from the name of this file:

cd "${0%/*}";

####### check cytosim

exe=./sim;

if [[ ! -x $exe ]]; then
	exe=../bin/sim;
fi

if [[ ! -x $exe ]]; then
	echo Error: missing executable $exe
	exit 1;
fi


####### make a new directory

n=0
while [[ $n -lt 9999 ]]; do
    dir=run$(printf "%04i" $n);
    if [[ ! -e $dir  ]]; then
        mkdir $dir;
        break;
    fi
    ((n+=1))
done

printf "\nThis will run in %s\n\n" $dir;

####### run cytosim

cp config.cym $dir/. || exit 1;
cd $dir;
../$exe;

####### salute

if [ $? == 0 ]; then
    printf "\nThank you for using Cytosim!\n\n";
else
    printf "\nCytosim did not terminate normally.\n";
    printf "Please report this problem to feedbackATcytosimDOTorg\n\n";
fi
