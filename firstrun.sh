#!/bin/bash
#PBS -N OCTA
#PBS -q q1b
#PBS -l ncpus=1
#PBS -j oe
#PBS -M s.matsu@phys-h.keio.ac.jp
#PMS -m b

ulimit -s unlimited

a="0.5373"
SCRIPT="make_octahedron_parameterfile01.sh"
cd $PBS_O_WORKDIR

for MM in 0.01 0.03 0.05 0.1 0.3 0.5 1.0; do 
	cd MM${MM}
	./${SCRIPT} 1 3 ${a} ${MM} 10 0.001 10 000 000
	sleep 1
	./gsugino01.exe > log &
	cd ../
done
wait
for MM in 0.01 0.03 0.05 0.1 0.3 0.5 1.0; do 
	cd MM${MM}
	./${SCRIPT} 0 3 ${a} ${MM} 10 0.01 90 000 000
	sleep 1
	./gsugino01.exe > log &
	cd ../
done
wait

for MM in 0.01 0.03 0.05 0.1 0.3 0.5 1.0; do 
	cd MM${MM}
	if [ "${MM}" == "0.01" ]; then
		./${SCRIPT} 0 3 ${a} ${MM} 10 0.01 1000 000 001
	else
		./${SCRIPT} 0 3 ${a} ${MM} 10 0.03 1000 000 001
	fi
	sleep 0.5
	./gsugino01.exe > log &
	cd ../
done

wait

