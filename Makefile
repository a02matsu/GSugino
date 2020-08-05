VER=04
#VER=debug
VER_CALCOBS=04
#FC=ifort
FC=mpiifort
#PARA=-DPARALLEL -DPARATEST
PARA=-DPARALLEL
PARA2=-DPARALLEL -DCOUNT_TIME
#PARA=-DNOPARALLEL
#FLAGS_IFORT=-mkl -fpp $(PARA) -CB -traceback -g 
FLAGS_CLUSTER=-mkl=cluster -fpp $(PARA) -CB -traceback -g 
#FLAGS_IFORT=-mkl -parallel -ipo
#FLAGS_IFORT=-mkl -fpp $(PARA) -O3 -ipo
#FLAGS_GCC=-llapack -lblas
# コンパイルのために順番が大事。下層ほど先に書く。 
SRCS=\
	SUN_generators.f90 \
	global_subroutines.f90 \
	matrix_functions.f90 \
	simplicial_complex.f90 \
	mt95.f90 \
	rational_algorithm.f90 \
	global_parameters.f90 \
	Dirac_operator.f90 \
	differential_Dirac.f90 \
	check_routines.f90 \
	simulation.f90 \
	parallel.f90 
OBJS=$(SRCS:.f90=.o)
#########################
SRC_MAIN=GSugino.f90  
OBJ_MAIN=GSugino.o
PROG=gsugino$(VER).exe
LIB=Pfapack_m02/libpfapack.a
########################
SRC_CALCOBSMAIN=calcobs.f90 
OBJ_CALCOBSMAIN=calcobs.o
PROG_CALCOBS=calcobs.exe
##
SRC_CALCOBSCOMM=initialization_calcobs.f90
OBJ_CALCOBSCOMM=$(SRC_CALCOBSCOMM:.f90=.o)
MEASUREMENT=Measurement
#########################
SRC_Dirac=writeDirac.f90  
OBJ_Dirac=writeDirac.o
PROG_Dirac=writeDirac.exe
#########################
SRC_Dinv2=calcDinv_CG.f90  
OBJ_Dinv2=calcDinv_CG.o
PROG_Dinv2=calcDinv_CG.exe
#########################
SRC_Dinv=calcDinv_PBLAS.f90  
OBJ_Dinv=calcDinv_PBLAS.o
PROG_Dinv=calcDinv.exe
#########################
SRC_TRPHI2=calc_trphi2.f90  
OBJ_TRPHI2=calc_trphi2.o
PROG_TRPHI2=calc_trphi2.exe
#########################
SRC_TRF2=calc_trf2.f90  
OBJ_TRF2=calc_trf2.o
PROG_TRF2=calc_trf2.exe
#########################
SRC_divJ_U1V=calc_divJ_U1V.f90  
OBJ_divJ_U1V=calc_divJ_U1V.o
PROG_divJ_U1V=calc_divJ_U1V.exe
#########################
SRC_localops=calc_local_operators.f90  
OBJ_localops=calc_local_operators.o
PROG_localops=calc_local_operators.exe
#########################
SRC_U1R=calc_exact_U1R.f90  
OBJ_U1R=calc_exact_U1R.o
PROG_U1R=calc_exact_U1R.exe
#########################
SRC_WriteConf=writeconfig.f90  
OBJ_WriteConf=writeconfig.o
PROG_WriteConf=writeconfig.exe
#########################
DIR_OBSERVABLES=Observables


#.SUFFIXES : .o .f90 # .oを作るときは必ず.f90から作るよ
.SUFFIXES : .f90 # .oを作るときは必ず.f90から作るよ
 
all:$(PROG) 

$(PROG): $(OBJS) $(OBJ_MAIN)
	$(FC) $(FLAGS_CLUSTER) -o $@  $(OBJS) $(OBJ_MAIN) $(LIB)
#########################################
obs:$(PROG_CALCOBS)
$(PROG_CALCOBS): $(OBJ_CALCOBSMAIN) $(OBJ_CALCOBSCOMM) $(OBJS)
	 $(FC) $(FLAGS_CLUSTER) -o $@ $(OBJ_CALCOBSCOMM) $(OBJS) $(OBJ_CALCOBSMAIN) $(LIB)
#########################################
U1V:$(PROG_divJ_U1V)

$(OBJ_divJ_U1V): $(MEASUREMENT)/FermionCorrelation_from_Dinv.f90 
$(PROG_divJ_U1V): $(OBJ_divJ_U1V) $(OBJS) $(OBJ_CALCOBSMAIN) 
	 $(FC) $(FLAGS_CLUSTER) -o $@ $(OBJ_CALCOBSCOMM) $(OBJS) $(OBJ_divJ_U1V) $(LIB)


#########################################
lops:$(PROG_localops)

$(OBJ_localops): $(MEASUREMENT)/FermionCorrelation_from_Dinv.f90 simulation.o $(OBJ_CALCOBSCOMM)
$(PROG_localops): $(OBJ_localops) $(OBJS) $(OBJ_CALCOBSMAIN) 
	 $(FC) $(FLAGS_CLUSTER) -o $@ $(OBJ_CALCOBSCOMM) $(OBJS) $(OBJ_localops) $(LIB)


writeconf:$(PROG_WriteConf)

$(PROG_WriteConf): $(OBJ_WriteConf)  $(OBJS) $(OBJ_CALCOBSCOMM) 
	$(FC) $(FLAGS_CLUSTER) -o $@ $(OBJ_CALCOBSCOMM) $(OBJS) $(OBJ_WriteConf) $(LIB)



#########################################
# moduleをコンパイルするときの依存性を解消
#structure_constant.o: structure_constant.f90
#	#$(FC) -c $<
#structure_constant.mod: structure_constant.f90 structure_constant.o
#	#@:
%.o: %.f90
	$(FC) $(FLAGS_CLUSTER) -c $<
%.mod: %.f90 %.o
	@true

# moduleの依存性
GSugino.o: \
  mt95.o \
  global_parameters.o \
  SUN_generators.o \
  simulation.o \
  parallel.o \
  check_routines.o 
global_parameters.o: \
  parallel.o \
  simplicial_complex.o \
  SUN_generators.o \
  set_theory_parameters.f90 \
  set_simulation_parameters.f90 \
  set_mpi_distribution.f90 \
  set_local_data.f90 \
  set_NZF.f90 \
  set_global_simplicial_complex.f90 \
  set_local_data.f90 \
  set_Remez_data.f90
global_subroutines.o: \
  global_parameters.o \
  parallel.o \
  mt95.o \
  SUN_generators.o \
  matrix_functions.o
initialization_calcobs.o: \
  global_parameters.f90
Dirac_operator.o: \
  global_parameters.o \
  global_subroutines.o \
  parallel.o \
  SUN_generators.o \
  matrix_functions.o \
  Dirac/prod_Dirac_site.f90 \
  Dirac/prod_Dirac_link1.f90 \
  Dirac/prod_Dirac_link2.f90 \
  Dirac/prod_Dirac_face1.f90 \
  Dirac/prod_Dirac_Omega.f90 \
  Dirac/prod_Dirac_Omega_m0.f90 \
  Dirac/prod_Dirac_Omega_adm.f90 \
  Dirac/prod_Dirac_mass.f90
differential_Dirac.o: \
  global_parameters.o \
  global_subroutines.o \
  parallel.o \
  SUN_generators.o \
  matrix_functions.o
simulation.o: \
  global_parameters.o \
  global_subroutines.o \
  mt95.o \
  parallel.o \
  Dirac_operator.o \
  differential_Dirac.o \
  matrix_functions.o \
  SUN_generators.o \
  rational_algorithm.o \
  MonteCarloSteps.f90 \
  hamiltonian.f90 \
  output.f90 \
  forces.f90 \
  check_QS.f90 \
  $(DIR_OBSERVABLES)/bosonic_action.f90 \
  $(DIR_OBSERVABLES)/fermionic_operators.f90 \
  $(DIR_OBSERVABLES)/div_rot.f90 \
  $(DIR_OBSERVABLES)/U1V_current.f90 \
  $(DIR_OBSERVABLES)/fermion_action_site.f90 \
  $(DIR_OBSERVABLES)/fermion_action_link1.f90 \
  $(DIR_OBSERVABLES)/fermion_action_link2.f90 \
  $(DIR_OBSERVABLES)/fermion_action_face1.f90 \
  $(DIR_OBSERVABLES)/fermion_action_face2.f90 \
  $(DIR_OBSERVABLES)/fermionic_face_lagrangian.f90 \
  $(DIR_OBSERVABLES)/compensators.f90 \
  $(DIR_OBSERVABLES)/phibar_compensator.f90 \
  $(DIR_OBSERVABLES)/Yphi_compensator.f90 \
  $(DIR_OBSERVABLES)/Yphibar_compensator.f90 \
  $(DIR_OBSERVABLES)/face_compensator.f90 \
  $(DIR_OBSERVABLES)/phichi.f90 \
  $(DIR_OBSERVABLES)/checkFF.f90 \
  $(DIR_OBSERVABLES)/trphi2.f90 \
  $(DIR_OBSERVABLES)/trf2.f90 \
  $(DIR_OBSERVABLES)/trivialWT.f90 \
  $(DIR_OBSERVABLES)/eigenvalues_of_Dirac.f90 \
  $(DIR_OBSERVABLES)/exact_U1R.f90 \
  $(DIR_OBSERVABLES)/mass_reweighting.f90 \
  $(DIR_OBSERVABLES)/WT_mass_contribution_site.f90 \
  $(DIR_OBSERVABLES)/WT_mass_contribution_link.f90 \
  $(DIR_OBSERVABLES)/WT_mass_contribution_face.f90 \
  $(DIR_OBSERVABLES)/make_Xi.f90 \
  $(DIR_OBSERVABLES)/Qfermion.f90 \
  $(DIR_OBSERVABLES)/QS_site.f90 \
  $(DIR_OBSERVABLES)/QS_link.f90 \
  $(DIR_OBSERVABLES)/QS_face.f90 \
  $(DIR_OBSERVABLES)/QS_3fermion_link.f90 \
  $(DIR_OBSERVABLES)/Xisite_Dinv.f90 \
  $(DIR_OBSERVABLES)/Xilink_Dinv.f90 \
  $(DIR_OBSERVABLES)/Xiface_Dinv.f90 
$(OBJ_CALCOBSMAIN): \
  global_parameters.o \
  differential_Dirac.o \
  simulation.o \
  initialization_calcobs.o \
  $(MEASUREMENT)/FermionCorrelation_from_Dinv.f90 \
  $(MEASUREMENT)/construct_Dirac.f90 \
  $(MEASUREMENT)/writeout_Dirac.f90 
$(OBJ_WriteConf): \
  global_parameters.o \
  simulation.o \
  initialization_calcobs.o \
  parallel.o 
(PROG_Dinv): \
  $(SRC_Dinv)


.PHONY: clean
clean:
	mv $(PROG) $(PROG).bak; rm -f *.o *.mod core 
