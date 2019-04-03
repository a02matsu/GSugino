#!/bin/bash

INFILE=$1
OUTFILE=$(echo $1 | sed "s/Dirac/Dinv/")

NMAT=3
num_sites=32
num_links=56
num_faces=26

NPROW=2
NPCOL=2
MB=2
NB=2

np=4

#echo $INFILE $NMAT
mpirun -n ${np} ./calcDinv.exe ${INFILE} ${OUTFILE} ${NMAT} ${num_sites} ${num_links} ${num_faces} ${NPROW} ${NPCOL} ${MB} ${NB} 


