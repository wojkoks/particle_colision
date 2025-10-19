#include <stdlib.h>
#include <stdio.h>
#include "particles.h"

void init_particle_array(ParticleArray *pa, int initial_capacity) {
    pa->count = initial_capacity;
    pa->capacity = initial_capacity;

    pa->array = malloc(sizeof(Particle) * pa->capacity);
    if (!pa->array) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(EXIT_FAILURE);
    }

    for (int i = 0; i < pa->count; i++) {
        pa->array[i].x = rand() % WIDTH;
        pa->array[i].y = rand() % HEIGHT;
        pa->array[i].vx = ((float)rand() / RAND_MAX - 0.5f) * 2 * V_MAX;
        pa->array[i].vy = ((float)rand() / RAND_MAX - 0.5f) * 2 * V_MAX;
        pa->array[i].color = rand() % 255;
    }
}

void free_particle_array(ParticleArray *pa) {
    if(pa->array) free(pa->array);
    pa->array = NULL;
    pa->count = 0;
    pa->capacity = 0;
}