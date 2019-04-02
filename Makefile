VER=01
#VER=debug
VER_CALCOBS=02
#FC=gfortran
#FC=ifort
FC=mpiifort
#PARA=-DPARALLEL -DPARATEST
PARA=-DPARALLEL
PARA2=-DPARALLEL -DCOUNT_TIME
#PARA=-DNOPARALLEL
#FLAGS_IFORT=-mkl -fpp $(PARA) -CB -traceback -g 
#FLAGS_IFORT=-mkl -parallel -ipo
FLAGS_IFORT=-mkl -fpp $(PARA) -O3 -ipo
FLAGS_GCC=-llapack -lblas
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
#########################
SRC_MAIN=GSugino.f90  
OBJ_MAIN=GSugino.o
PROG=gsugino$(VER).exe
LIB=libpfapack.a
#########################
SRC_OBS=calcobs.f90  
OBJ_OBS=calcobs.o
PROG_OBS=calcobs.exe
#########################
SRC_Dirac=writeDirac.f90  
OBJ_Dirac=writeDirac.o
PROG_Dirac=writeDirac.exe
#########################
SRC_Dinv2=writeDinv.f90  
OBJ_Dinv2=writeDinv.o
PROG_Dinv2=writeDinv.exe
#########################
SRC_Dinv=calcDinv.f90  
OBJ_Dinv=calcDinv.o
PROG_Dinv=calcDinv.exe
#########################
LIB=libpfapack.a


#.SUFFIXES : .o .f90 # .oを作るときは必ず.f90から作るよ
.SUFFIXES : .f90 # .oを作るときは必ず.f90から作るよ
 
all:$(PROG) 

$(PROG): $(OBJS) $(OBJ_MAIN)
ifeq ($(FC),gfortran)
	$(FC) -O2 $(FLAGS_GCC) -o $@ $(OBJS) $(OBJ_MAIN) $(LIB)
else
	$(FC) $(FLAGS_IFORT) -o $@ $(OBJS) $(OBJ_MAIN) $(LIB)
endif

obs:$(PROG_OBS)

$(PROG_OBS): $(OBJ_OBS) $(OBJ_MAIN)
ifeq ($(FC),gfortran)
	$(FC) -O2 $(FLAGS_GCC) -o $@ $(OBJS) $(OBJ_OBS) $(LIB)
else
	$(FC) $(FLAGS_IFORT) -o $@ $(OBJS) $(OBJ_OBS) $(LIB)
endif

dirac:$(PROG_Dirac)

$(PROG_Dirac): $(OBJ_Dirac) $(OBJ_MAIN)
ifeq ($(FC),gfortran)
	$(FC) -O2 $(FLAGS_GCC) -o $@ $(OBJS) $(OBJ_Dirac) $(LIB)
else
	$(FC) $(FLAGS_IFORT) -o $@ $(OBJS) $(OBJ_Dirac) $(LIB)
endif

dinv:$(PROG_Dinv)

$(PROG_Dinv): $(SRC_Dinv)
	mpiifort -mkl=cluster -CB -traceback -g $(SRC_Dinv) -o $(PROG_Dinv)

dinv2:$(PROG_Dinv2)

$(PROG_Dinv2): $(OBJ_Dinv2) $(OBJ_MAIN)
ifeq ($(FC),gfortran)
	$(FC) -O2 $(FLAGS_GCC) -o $@ $(OBJS) $(OBJ_Dinv2) $(LIB)
else
	$(FC) $(FLAGS_IFORT) -o $@ $(OBJS) $(OBJ_Dinv2) $(LIB)
endif

# moduleをコンパイルするときの依存性を解消
#structure_constant.o: structure_constant.f90
#	#$(FC) -c $<
#structure_constant.mod: structure_constant.f90 structure_constant.o
#	#@:
%.o: %.f90
ifeq ($(FC),gfortran)
	$(FC) -c $<
else
	$(FC) $(FLAGS_IFORT) -c $<
endif
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
  bosonic_action.f90 \
  WT_identities.f90 \
  divK3.f90 \
  divK4.f90 \
  rotK1.f90 \
  rotK2.f90 \
  U1V_current.f90
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


.PHONY: clean
clean:
	mv $(PROG) $(PROG).bak; rm -f *.o *.mod core 
