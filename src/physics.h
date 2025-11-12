#pragma once
#include <SDL2/SDL.h>

typedef struct
{
    float x, y;
} Vec2;

typedef struct
{
    Vec2 pos;
    Vec2 vel;
    float radius;
    float mass;
    Uint8 r, g, b;
} Particle;

typedef struct
{
    int width, height;
    int count;
    Particle *p;
} World;

void world_init(World *w, int width, int height, int count, unsigned seed,
                int rmin, int rmax, float max_speed);
void world_free(World *w);
void world_reset(World *w, unsigned seed);
void world_step(World *w, float dt);