#ifndef PARTICLES_H
#define PARTICLES_H

#define WIDTH  10
#define HEIGHT 10
#define V_MAX  1.0f
#define DT     0.01f

typedef struct {
    float x, y;
    float vx, vy;
    int color;
} Particle;


typedef struct {
    Particle *array;
    int count;
    int capacity;
} ParticleArray;

void init_particle_array(ParticleArray *pa, int initial_capacity);
void free_particle_array(ParticleArray *pa);

#endif