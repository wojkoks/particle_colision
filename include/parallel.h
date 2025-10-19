#ifndef PARALLEL_H
#define PARALLEL_H

#include "particles.h"
#include "grid.h"

void update_positions_openmp(ParticleArray *pa);
void update_positions_mpi(ParticleArray *pa);
void update_positions_gpu(ParticleArray *pa);

#endif