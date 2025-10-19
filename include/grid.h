#ifndef GRID_H
#define GRID_H

#include "particles.h"

typedef struct {
    int cols, rows;
    ParticleArray **cells;
} Grid;

Grid *create_grid(int cols, int rows);
void free_grid(Grid *grid);
void assign_particles_to_grid(Grid *grid, ParticleArray *pa);
void check_collisions_in_grid(Grid *grid);

#endif