VER=01
VER_CALCOBS=02
#FC=gfortran
FC=ifort
#FLAGS_IFORT=-mkl -CB -traceback -g 
#FLAGS_IFORT=-mkl -parallel -ipo
FLAGS_IFORT=-mkl -O2 
FLAGS_GCC=-llapack -lblas
# コンパイルのために順番が大事。下層ほど先に書く。 
SRCS=SUN_generators.f90 \
	 matrix_functions.f90 \
	 simplicial_complex.f90 \
	 mt95.f90 \
	 rational_algorithm.f90 \
	 global_parameters.f90 \
	 global_subroutines.f90 \
	 initialization.f90 \
	 Dirac_operator.f90 \
	 differential_Dirac.f90 \
	 hamiltonian.f90 \
	 forces.f90 \
	 observables.f90 \
	 output.f90 \
	 simulation.f90 
OBJS=$(SRCS:.f90=.o)
#########################
SRC_MAIN=GSugino.f90  
OBJ_MAIN=GSugino.o
PROG=gsugino$(VER).exe
#########################
SRC_CALCOBS=calcobs.f90
OBJ_CALCOBS=calcobs.o
PROG_CALCOBS=calcobs$(VER_CALCOBS).exe
#########################
SRC_CALCCOMP=calc_compensator.f90
OBJ_CALCCOMP=calc_compensator.o
PROG_CALCCOMP=calc_compensator$(VER_CALCOBS).exe
#########################
SRC_CALCEIGEN=calceigen.f90
OBJ_CALCEIGEN=calceigen.o
PROG_CALCEIGEN=calceigen$(VER_CALCOBS).exe
#########################
SRC_CALCWT=calcWTIZ.f90
OBJ_CALCWT=calcWTIZ.o
PROG_CALCWT=calcWTIZ$(VER_CALCOBS).exe
#########################
SRC_WRITE_DOWN_CONFIG=write_down_config.f90
OBJ_WRITE_DOWN_CONFIG=write_down_config.o
PROG_WRITE_DOWN_CONFIG=write_down_config.exe
########################
LIB=libpfapack.a


#.SUFFIXES : .o .f90 # .oを作るときは必ず.f90から作るよ
.SUFFIXES : .f90 # .oを作るときは必ず.f90から作るよ
 
all:$(PROG) $(PROG_CALCOBS) $(PROG_CALCCOMP) $(PROG_CALCEIGEN) $(PROG_WRITE_DOWN_CONFIG) $(PROG_CALCWT)

#main:$(PROG)

#calcobs:$(PROG_CALCOBS)

#write_down_config:$(PROC_WRITE_DOWN_CONFIG)

$(PROG): $(OBJS) $(OBJ_MAIN)
ifeq ($(FC),gfortran)
	$(FC) -O2 $(FLAGS_GCC) -o $@ $(OBJS) $(OBJ_MAIN) $(LIB)
else
	$(FC) -O2 $(FLAGS_IFORT) -o $@ $(OBJS) $(OBJ_MAIN) $(LIB)
endif

$(PROG_CALCOBS): $(OBJ_CALCOBS) $(OBJS) 
ifeq ($(FC),gfortran)
	$(FC) -O2 $(FLAGS_GCC) -o $@ $(OBJS) $(OBJ_CALCOBS) $(LIB)
else
	$(FC) -O2 $(FLAGS_IFORT) -o $@ $(OBJS) $(OBJ_CALCOBS) $(LIB)
endif


$(PROG_CALCCOMP): $(OBJ_CALCCOMP) $(OBJS) 
ifeq ($(FC),gfortran)
	$(FC) -O2 $(FLAGS_GCC) -o $@ $(OBJS) $(OBJ_CALCCOMP) $(LIB)
else
	$(FC) -O2 $(FLAGS_IFORT) -o $@ $(OBJS) $(OBJ_CALCCOMP) $(LIB)
endif

$(PROG_CALCEIGEN): $(OBJ_CALCEIGEN) $(OBJS) 
ifeq ($(FC),gfortran)
	$(FC) -O2 $(FLAGS_GCC) -o $@ $(OBJS) $(OBJ_CALCEIGEN) $(LIB)
else
	$(FC) -O2 $(FLAGS_IFORT) -o $@ $(OBJS) $(OBJ_CALCEIGEN) $(LIB)
endif

$(PROG_CALCWT) : $(OBJ_CALCWT) $(OBJS)
ifeq ($(FC),gfortran)
	$(FC) $(FLAGS_GCC) -o $@ $(OBJS) $(OBJ_CALCWT) $(LIB)
else
	$(FC) -O2 $(FLAGS_IFORT) -o $@ $(OBJS) $(OBJ_CALCWT) $(LIB)
endif

$(PROG_WRITE_DOWN_CONFIG) : $(OBJ_WRITE_DOWN_CONFIG) $(OBJS)
ifeq ($(FC),gfortran)
	$(FC) $(FLAGS_GCC) -o $@ $(OBJS) $(OBJ_WRITE_DOWN_CONFIG) $(LIB)
else
	$(FC) -O2 $(FLAGS_IFORT) -o $@ $(OBJS) $(OBJ_WRITE_DOWN_CONFIG) $(LIB)
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
  initialization.o \
  simulation.o \
  simplicial_complex.o \
  matrix_functions.o \
  global_subroutines.o \
  SUN_generators.o 
global_parameters.o: \
  simplicial_complex.o \
  SUN_generators.o
global_subroutines.o: \
  global_parameters.o \
  mt95.o \
  SUN_generators.o \
  matrix_functions.o
initialization.o: \
  global_parameters.o \
  mt95.o \
  SUN_generators.o \
  matrix_functions.o \
  global_subroutines.o
Dirac_operator.o: \
  global_parameters.o \
  global_subroutines.o \
  SUN_generators.o \
  matrix_functions.o
differential_Dirac.o: \
  global_parameters.o \
  global_subroutines.o \
  SUN_generators.o \
  matrix_functions.o
hamiltonian.o: \
  global_parameters.o \
  global_subroutines.o \
  SUN_generators.o \
  matrix_functions.o \
  rational_algorithm.o \
  Dirac_operator.o
forces.o: \
  global_parameters.o \
  global_subroutines.o \
  SUN_generators.o \
  matrix_functions.o \
  rational_algorithm.o \
  Dirac_operator.o \
  differential_Dirac.o
observables.o: \
  global_parameters.o \
  Dirac_operator.o \
  SUN_generators.o \
  matrix_functions.o
output.o: \
  global_parameters.o \
  observables.o
simulation.o: \
  global_parameters.o \
  global_subroutines.o \
  hamiltonian.o \
  forces.o \
  mt95.o \
  Dirac_operator.o \
  output.o \
  matrix_functions.o \
  SUN_generators.o \
  rational_algorithm.o

$(OBJ_CALCOBS): $(OBJS)

$(OBJ_WRITE_DOWN_CONFIG): $(OBJS)


# moduleの依存性
#initialization.o: mt95.o matrix_functions

#%.o: %.mod
#$(OBJS):$(SRCS)
# メインコードの依存性
#$(MAIN_OBJ):$(OTHER_OJS) #.f90.o: #	$(FC) -c $< 

.PHONY: clean
clean:
	rm -f *.o *.mod $(PROG) core $(PROG_CALCOBS)
