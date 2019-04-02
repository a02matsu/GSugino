#!/bin/bash

INFILE=$1
OUTFILE=$(echo $1 | sed "s/Dirac/Dinv/")

NMAT=3
num_sites=8
num_links=12
num_faces=6

NPROW=2
NPCOL=2
MB=1
NB=1

np=4

#echo $INFILE $NMAT
mpirun -n ${np} ./calcDinv.exe ${INFILE} ${OUTFILE} ${NMAT} ${num_sites} ${num_links} ${num_faces} ${NPROW} ${NPCOL} ${MB} ${NB}


