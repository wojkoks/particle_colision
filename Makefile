CFLAGS  = -std=c11 -O2 -Wall -Wextra -Wpedantic -fopenmp $(shell pkg-config --cflags sdl2)
LDFLAGS = $(shell pkg-config --libs sdl2) -lm -fopenmp
MPICC   = mpicc
MPIFLAGS = -std=c11 -O2 -Wall -Wextra -Wpedantic -fopenmp $(shell pkg-config --cflags sdl2)
MPILDFLAGS = $(shell pkg-config --libs sdl2) -lm -fopenmp
NVCC    = nvcc
CUDA_ARCH ?= sm_86
CUDAFLAGS = -std=c++14 -O2 -arch=$(CUDA_ARCH) $(shell pkg-config --cflags sdl2)
CUDALDFLAGS = $(shell pkg-config --libs sdl2) -lm
NVCC_PATH := $(shell command -v $(NVCC) 2>/dev/null)
NVCC_CCBIN ?=
ifdef NVCC_CCBIN
NVCC_CCBIN_FLAG := -ccbin $(NVCC_CCBIN)
endif

ifeq ($(filter particles_cuda,$(MAKECMDGOALS)),particles_cuda)
ifeq ($(NVCC_PATH),)
$(error "nvcc not found. Install CUDA toolkit or set NVCC=/path/to/nvcc")
endif
endif

.PHONY: all clean

all: particles particles_mpi

# OpenMP version
particles: src/main.c src/physics.c src/draw.c
	$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS)

# MPI version
particles_mpi: src/main_mpi.c src/physics_mpi.c src/draw.c
	$(MPICC) $(MPIFLAGS) $^ -o $@ $(MPILDFLAGS)

# CUDA version (GPU)
particles_cuda: src/main.c src/physics.c src/physics_cuda.cu src/draw.c
	$(NVCC) $(CUDAFLAGS) $(NVCC_CCBIN_FLAG) -DUSE_CUDA -Xcompiler "-O2 -Wall -Wextra -Wpedantic" $^ -o $@ $(CUDALDFLAGS)

clean:
	rm -f particles particles_mpi particles_cuda

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
