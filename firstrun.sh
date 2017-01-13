$!/bin/bash
#PBS -N test_run
#PBS -q q1b
#PBS -l ncpus=1
#PBS -j oe

cd $PBS_O_WORKDIR

for MM in 0.01 0.03 0.05 0.1 0.3 0.5 1.0; do 
	cd MM${MM}
	./make_tetra_parameterfile01.sh 1 3 0.7598 ${MM} 10 0.001 10 000 000
	sleep 1
	./gsugino01.exe > log &
	cd ../
done

wait

