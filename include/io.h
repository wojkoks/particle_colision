#ifndef IO_H
#define IO_H

#include "particles.h"

void save_frame_ppm(ParticleArray *pa, int frame_number);
void print_fps(double fps);

#endif