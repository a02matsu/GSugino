CC = g++
CFLAGS = -O3
LFLAGS = ../../libpfapack.a -llapack -lblas -lgfortran

all: example_dense.out example_band.out

example_dense.out: example_dense.cc
	$(CC) $(CFLAGS) example_dense.cc -o example_dense.out $(LFLAGS)

example_band.out: example_band.cc
	$(CC) $(CFLAGS) example_band.cc -o example_band.out $(LFLAGS)

clean:
	rm -f *.out