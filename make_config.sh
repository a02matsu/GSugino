#!/bin/bash
#PBS -N TETRA
#PBS -q q1h
#PBS -l ncpus=1
#PBS -M s.matsu@phys-h.keio.ac.jp
#PMS -m b

cd $PBS_O_WORKDIR

for MM in 0.01 0.03 0.05 0.1 0.3 0.5 1.0; do 
	cd MM${MM}
	if [ "${MM}" == "0.01" ]; then
		./make_tetra_parameterfile01.sh 0 3 0.7598 ${MM} 10 0.01 1000 000 001
	else
		./make_tetra_parameterfile01.sh 0 3 0.7598 ${MM} 10 0.03 1000 000 001
	fi
	sleep 0.5
	./gsugino01.exe > log &
	cd ../
done

wait

