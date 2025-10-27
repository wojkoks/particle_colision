#pragma once
#include <SDL2/SDL.h>

typedef struct
{
    float x, y;
} Vec2;

typedef struct
{
    Vec2 pos;      // position in pixels
    Vec2 vel;      // velocity in px/s
    float radius;  // radius in px
    float mass;    // mass
    Uint8 r, g, b; // color
} Particle;

typedef struct
{
    int width, height;
    int count;   // current number of particles
    Particle *p; // array of particles (size >= count)
} World;

void world_init(World *w, int width, int height, int count, unsigned seed,
                int rmin, int rmax, float max_speed);
void world_free(World *w);
void world_reset(World *w, unsigned seed);
void world_step(World *w, float dt);