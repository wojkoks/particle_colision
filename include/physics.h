#ifndef PHYSICS_H
#define PHYSICS_H

#include "particles.h"

void update_positions(ParticleArray *pa);
void handle_collisions(ParticleArray *pa);
void bounce_from_walls(Particle *p);

#endif