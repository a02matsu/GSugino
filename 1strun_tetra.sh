#!/bin/bash

POL=$(echo $0 | cut -d _ -f 2 | cut -d . -f 1)
for MM in 0.01 0.03 0.05 0.1 0.3 0.5 1.0; do 
./run_${POL}.sh ${MM} 0.0001 10 000 000 10 1

done
