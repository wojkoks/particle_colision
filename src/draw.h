#pragma once
#include <SDL2/SDL.h>

void draw_filled_circle(SDL_Renderer *r, int cx, int cy, int radius);
void draw_grid(SDL_Renderer *r, int width, int height, int cols, int rows);