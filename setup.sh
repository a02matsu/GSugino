#!/bin/bash


PARASCRIPT="make_"$1"_parameterfile01.sh"

cp ~/GitHub/GSugino/make_tetra_parameterfile01.sh ./${PARASCRIPT}
for MM in 0.01 0.03 0.05 0.1 0.3 0.5 1.0; do 
	mkdir MM${MM}
	cd MM${MM}
	ln -s ../${PARASCRIPT} .
	ln -s ~/GitHub/GSugino/gsugino01.exe .
	cp ~/GitHub/GSugino/remez_Q24* .
	cp ~/GitHub/GSugino/$1.dat .
	cd ..
done



