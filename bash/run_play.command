#! /bin/bash

####### go to directory extracted from the name of this file:

cd "${0%/*}";

####### find executable

exe=./play;

if [[ ! -x $exe ]]; then
	exe=./cytosim;
fi

if [[ ! -x $exe ]]; then
	exe=bin/play;
fi

if [[ ! -x $exe ]]; then
	echo Error: missing executable $exe
	exit 1;
fi

####### run

$exe;

