VER=01
VER_CALCOBS=02
#FC=gfortran
#FC=ifort
FC=mpiifort
#PARA=-DPARALLEL -DPARATEST
PARA=-DPARALLEL
#PARA=-DNOPARALLEL
FLAGS_IFORT=-mkl -fpp $(PARA) -CB -traceback -g 
#FLAGS_IFORT=-mkl -parallel -ipo
#FLAGS_IFORT=-mkl -O2 
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


#.SUFFIXES : .o .f90 # .oを作るときは必ず.f90から作るよ
.SUFFIXES : .f90 # .oを作るときは必ず.f90から作るよ
 
all:$(PROG) 

$(PROG): $(OBJS) $(OBJ_MAIN)
ifeq ($(FC),gfortran)
	$(FC) -O2 $(FLAGS_GCC) -o $@ $(OBJS) $(OBJ_MAIN) $(LIB)
else
	$(FC) -O2 $(FLAGS_IFORT) -o $@ $(OBJS) $(OBJ_MAIN) $(LIB)
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
  set_local_data.f90
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
#initialization.o: \
#  global_parameters.o \
#  mt95.o \
#  parallel.o \
#  SUN_generators.o \
#  matrix_functions.o \
#  global_subroutines.o
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
  forces.f90 
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
	rm -f *.o *.mod $(PROG) core $(PROG_CALCOBS)
