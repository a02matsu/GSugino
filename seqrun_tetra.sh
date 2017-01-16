#!/bin/bash

DTAU=$1
ITE=$2
OLD=$3
tmp=`expr $OLD + 1`
NEW=`printf "%03d\n" $tmp`
POL=$(echo $0 | cut -d _ -f 2 | cut -d . -f 1)
for MM in 0.01 0.03 0.05 0.1 0.3 0.5 1.0; do 
./run_${POL}.sh ${MM} ${DTAU} 10 $OLD $NEW $ITE

done
