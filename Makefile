CFLAGS  = -std=c11 -O2 -Wall -Wextra -Wpedantic -fopenmp $(shell pkg-config --cflags sdl2)
LDFLAGS = $(shell pkg-config --libs sdl2) -lm -fopenmp
MPICC   = mpicc
MPIFLAGS = -std=c11 -O2 -Wall -Wextra -Wpedantic -fopenmp $(shell pkg-config --cflags sdl2)
MPILDFLAGS = $(shell pkg-config --libs sdl2) -lm -fopenmp

.PHONY: all clean

all: particles particles_mpi

# OpenMP version
particles: src/main.c src/physics.c src/draw.c
	$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS)

# MPI version
particles_mpi: src/main_mpi.c src/physics_mpi.c src/draw.c
	$(MPICC) $(MPIFLAGS) $^ -o $@ $(MPILDFLAGS)

clean:
	rm -f particles particles_mpi

run-mpi-optimized:
	export OMP_NUM_THREADS=2; \
	export OMP_WAIT_POLICY=PASSIVE; \
	export OMP_DYNAMIC=FALSE; \
	mpirun -np 4 ./particles_mpi

run-mpi-single-thread:
	export OMP_NUM_THREADS=1; \
	mpirun -np 4 ./particles_mpi

run-omp:
	export OMP_NUM_THREADS=4; \
	./particles