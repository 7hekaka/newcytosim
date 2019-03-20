#! /bin/bash

####### go to directory extracted from the name of this file:

cd "${0%/*}";

####### make directory

n=0
while [[ $n -lt 9999 ]]; do
    dir=run$(printf "%04i" $n);
    if [[ ! -e $dir  ]]; then
        mkdir $dir;
        break;
    fi
    ((n+=1))
done

####### save run

mv *.cmo $dir;
cp config.cym $dir;

####### salute

printf "\nRun was saved in $dir\n\n";

