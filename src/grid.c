#include <stdlib.h>
#include "grid.h"

Grid *create_grid(int cols, int rows) {
    Grid *g = malloc(sizeof(Grid));
    g->cols = cols;
    g->rows = rows;
    g->cells = calloc(cols * rows, sizeof(ParticleArray*));
    return g;
}

void free_grid(Grid *grid) {
    free(grid->cells);
    free(grid);
}

void assign_particles_to_grid(Grid *grid, ParticleArray *pa) {
    (void)grid;
    (void)pa;
}

void check_collisions_in_grid(Grid *grid) {
    (void)grid;
}