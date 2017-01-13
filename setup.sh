#!/bin/bash

CODEDIR='~/GitHub/GSugino'
PARASCRIPT="make_"$1"_parameterfile01.sh"

mkdir $1
cd $1
cp ${CODEDIR}/make_tetra_parameterfile01.sh ./${PARASCRIPT}
cp ${CODEDIR}/firstrun.sh .
cp ${CODEDIR}/make_config.sh .
for MM in 0.01 0.03 0.05 0.1 0.3 0.5 1.0; do 
	mkdir MM${MM}
	cd MM${MM}
	ln -s ../${PARASCRIPT} .
	ln -s ${CODEDIR}/gsugino01.exe .
	cp ${CODEDIR}/remez_Q24* .
	cp ${CODEDIR}/$1.dat .
	cd ..
done



