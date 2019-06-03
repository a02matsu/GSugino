VER=01
#VER=debug
VER_CALCOBS=02
#FC=ifort
FC=mpiifort
#PARA=-DPARALLEL -DPARATEST
PARA=-DPARALLEL
PARA2=-DPARALLEL -DCOUNT_TIME
#PARA=-DNOPARALLEL
FLAGS_IFORT=-mkl -fpp $(PARA) -CB -traceback -g 
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
	check_routines.f90 \
	Dirac_operator.f90 \
	differential_Dirac.f90 \
	simulation.f90 \
	parallel.f90 
OBJS=$(SRCS:.f90=.o)
DIR_OBS=Observables
#########################
SRC_MAIN=GSugino.f90  
OBJ_MAIN=GSugino.o
PROG=gsugino$(VER).exe
LIB=Pfapack_m02/libpfapack.a
########################
SRC_OBSMAIN=calcobs.f90 
OBJ_OBSMAIN=calcobs.o
PROG_OBS=calcobs.exe
##
SRC_OBS=initialization_calcobs.f90
OBJ_OBS=$(SRC_OBS:.f90=.o)
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
PROG_Dinv=calcDinv_PBLAS.exe
#########################
SRC_CORR=calc_correlations.f90  
OBJ_CORR=calc_correlations.o
PROG_CORR=calc_correlations.exe
#########################


#.SUFFIXES : .o .f90 # .oを作るときは必ず.f90から作るよ
.SUFFIXES : .f90 # .oを作るときは必ず.f90から作るよ
 
all:$(PROG) 

$(PROG): $(OBJS) $(OBJ_MAIN)
	$(FC) $(FLAGS_IFORT) -o $@  $(OBJS) $(OBJ_MAIN) $(LIB)

obs:$(PROG_OBS)
#########################################
#$(PROG_OBS): $(OBJ_OBS) $(OBJ_MAIN)
$(PROG_OBS): $(OBJ_OBS) $(OBJ_OBSMAIN)
	 $(FC) $(FLAGS_IFORT) -o $@ $(OBJ_OBS) $(OBJS) $(OBJ_OBSMAIN) $(LIB)
#########################################
dirac:$(PROG_Dirac)

$(PROG_Dirac): $(OBJ_Dirac) $(OBJ_MAIN)
	$(FC) $(FLAGS_IFORT) -o $@ $(OBJS) $(OBJ_Dirac) $(LIB)
#########################################
dinv:$(PROG_Dinv)

$(PROG_Dinv): $(SRC_Dinv)
	$(FC) $(FLAGS_IFORT) -o $@ $(SRC_Dinv) 
#########################################
dinv2:$(PROG_Dinv2)

$(PROG_Dinv2): $(OBJ_Dinv2) $(OBJ_MAIN)
	$(FC) $(FLAGS_IFORT) -o $@ $(OBJS) $(OBJ_Dinv2) $(LIB)
#########################################
corr:$(PROG_CORR)

$(PROG_CORR): $(OBJ_CORR) $(OBJ_OBS) $(OBJS) $(OBJ_CORR)
	$(FC) $(FLAGS_IFORT) -o $@ $(OBJS) $(OBJ_OBS) $(OBJ_CORR) $(LIB)
#########################################
# moduleをコンパイルするときの依存性を解消
#structure_constant.o: structure_constant.f90
#	#$(FC) -c $<
#structure_constant.mod: structure_constant.f90 structure_constant.o
#	#@:
%.o: %.f90
	$(FC) $(FLAGS_IFORT) -c $<
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
check_routines.o: \
  global_parameters.o \
  global_subroutines.o \
  parallel.o \
  matrix_functions.o \
  mt95.o \
  simplicial_complex.o
simulation.o: \
  global_parameters.o \
  global_subroutines.o \
  mt95.o \
  parallel.o \
  Dirac_operator.o \
  matrix_functions.o \
  SUN_generators.o \
  rational_algorithm.o \
  MonteCarloSteps.f90 \
  hamiltonian.f90 \
  observables.f90 \
  output.f90 \
  forces.f90 \
  $(DIR_OBS)/bosonic_action.f90 \
  $(DIR_OBS)/fermionic_operators.f90 \
  $(DIR_OBS)/WT_identities.f90 \
  $(DIR_OBS)/divK3.f90 \
  $(DIR_OBS)/divK4.f90 \
  $(DIR_OBS)/rotK1.f90 \
  $(DIR_OBS)/rotK2.f90 \
  $(DIR_OBS)/U1V_current.f90 \
  $(DIR_OBS)/compensators.f90 \
  $(DIR_OBS)/phichi.f90 \
  $(DIR_OBS)/checkFF.f90 \
  $(DIR_OBS)/local_operators.f90
Dirac_operator.o: \
  global_parameters.o \
  global_subroutines.o \
  parallel.o \
  SUN_generators.o \
  matrix_functions.o
differential_Dirac.o: \
  global_parameters.o \
  global_subroutines.o \
  parallel.o \
  SUN_generators.o \
  matrix_functions.o
initialization_calcobs.o: \
  global_parameters.f90



.PHONY: clean
clean:
	mv $(PROG) $(PROG).bak; rm -f *.o *.mod core 
