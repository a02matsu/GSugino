#!/bin/bash

# usage: 
# ./make_paramterfile**.sh NewConf(0/1) LatticeSpacing mass^2 Ntau Dtau NumIte OldNum

PARAFILE="parameters.dat"

SimpFile="tetra.dat"
ABFile="ab_tetra.dat"

NewConf=$1
NMAT=$2
LatticeSpacing=$3
MassSq=$4
Ntau=$5
Dtau=$6
NumIte=$7
OLDNUM=$8
NEWNUM=$9

#PhysMass2=$( echo "scale=6; ${MassSq}*${LatticeSpacing}*${LatticeSpacing}" | bc)

if [ -z "$NUWNUM" ]; then 
	if [ "$NewConf" = "0" ] ; then
  	tmp=`expr $OLDNUM + 1`
  	NEWNUM=`printf "%03d\n" $tmp`
	else
  	NEWNUM=${OLDNUM}
	fi
fi

echo "${NMAT}			! NMAT" > ${PARAFILE}
echo "${LatticeSpacing}d0			! lattice spacing" >> ${PARAFILE}
echo "${SimpFile}		! data of the simplicial complex" >> ${PARAFILE}
echo "${ABFile}	! data of alpha and beta" >> ${PARAFILE}
echo "0				! test_mode ; 0:Simulation mode 1:Test mode, " >> ${PARAFILE}
echo "${NewConf}				! new_config ; 0:Old Config 1:New Config" >> ${PARAFILE}
echo "0				! fix_seed ; 0:previous value 1:fix 2:system time" >> ${PARAFILE}
echo "0				! reset_ite ; 0:continue total ite, 1:reset total ite" >> ${PARAFILE}
echo "0				! read_alpha; 0:alpha,beta=1d0 1:read alpha/beta from file" >> ${PARAFILE}
echo "1				! save_med_step; step to write out configuration" >> ${PARAFILE}
echo "100			! save_config_step; writedown configfile by this step" >> ${PARAFILE}
echo "1				! obs_step; step to compute observables" >> ${PARAFILE}
echo "42342			! seed (works when fix_seed=1)" >> ${PARAFILE}
echo "1				! m_omega ;integer to avoid vacuum degeneracy (mm>=NMAT/4)" >> ${PARAFILE}
echo "${MassSq}d0 		! mass_phi2 : physical mass square of \Phi" >> ${PARAFILE}
echo "0.0d0       	! mass_f; fermion mass" >> ${PARAFILE}
echo "1d4		        ! Remez_factor4 ; range of remez approx( min*factor .. max*factor) " >> ${PARAFILE}
echo "1d4		        ! Remez_factor8 ; range of remez approx( min*factor .. max*factor) " >> ${PARAFILE}
echo "1d-12	    	! epsilon ; range to stop CG solver" >> ${PARAFILE}
echo "10000	    	! CG_max ; maximum number of CG iteration" >> ${PARAFILE}
echo "${NumIte}     	! num_ite ; number of iterations" >> ${PARAFILE}
echo "${Ntau}		! NTAU" >> ${PARAFILE}
echo "${Dtau}d0			! DTAU base" >> ${PARAFILE}
echo "6d0				! R_phi ; Dtau_phi=R_phi * Dtau" >> ${PARAFILE}
echo "1d0				! R_A   ; Dtau_A=R_A * Dtau" >> ${PARAFILE}
echo "conSU${NMAT}m1a${LatticeSpacing}PM${MassSq}_${OLDNUM}.dat ! input configuration " >> ${PARAFILE}
echo "outSU${NMAT}m1a${LatticeSpacing}PM${MassSq}_${NEWNUM}.dat ! output data" >> ${PARAFILE}
echo "conSU${NMAT}m1a${LatticeSpacing}PM${MassSq}_${NEWNUM}.dat ! output configuration " >> ${PARAFILE}
echo "medSU${NMAT}m1a${LatticeSpacing}PM${MassSq}_${NEWNUM}.dat ! intermediate configuration " >> ${PARAFILE}

