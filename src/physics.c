#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <SDL2/SDL.h>
#include "physics.h"

#ifdef _OPENMP
#include <omp.h>
#endif

void world_init(World *w, int width, int height, int count, unsigned seed, int rmin, int rmax, float max_speed)
{
    w->width = width;
    w->height = height;
    w->count = count;
    w->p = (Particle *)malloc(count * sizeof(Particle));
    srand(seed);

    for (int i = 0; i < count; i++)
    {
        w->p[i].radius = rand() % (rmax - rmin + 1) + rmin;
        w->p[i].pos.x = w->p[i].radius + (float)(rand()) / RAND_MAX * (width - 2 * w->p[i].radius);
        w->p[i].pos.y = w->p[i].radius + (float)(rand()) / RAND_MAX * (height - 2 * w->p[i].radius);
        w->p[i].vel.x = ((float)rand() / RAND_MAX) * max_speed * (rand() % 2 == 0 ? 1 : -1);
        w->p[i].vel.y = ((float)rand() / RAND_MAX) * max_speed * (rand() % 2 == 0 ? 1 : -1);
        w->p[i].mass = 3.14f * w->p[i].radius * w->p[i].radius;
        w->p[i].r = rand() % 256;
        w->p[i].g = rand() % 256;
        w->p[i].b = rand() % 256;
    }
    
    // Debug: calculate velocity stats
    float vel_sum = 0, vel_max = 0, vel_min = 999;
    for (int i = 0; i < count; i++)
    {
        float speed = sqrtf(w->p[i].vel.x * w->p[i].vel.x + w->p[i].vel.y * w->p[i].vel.y);
        vel_sum += speed;
        if (speed > vel_max) vel_max = speed;
        if (speed < vel_min) vel_min = speed;
    }
    printf("[OMP] Velocity stats: avg=%.2f min=%.2f max=%.2f (expected_max=%.2f)\n", 
           vel_sum / count, vel_min, vel_max, max_speed * 1.414f);
}

void world_free(World *w)
{
    free(w->p);
    w->p = NULL;
}

void world_reset(World *w, unsigned seed)
{
    int width = w->width;
    int height = w->height;
    int count = w->count;
    float max_speed = 120.0f;
    int rmin = 3;
    int rmax = 8;

    // Free current particles
    world_free(w);

    // Reinitialize with same parameters but new seed
    world_init(w, width, height, count, seed, rmin, rmax, max_speed);
}

void world_step(World *w, float dt)
{
    // === OPENMP: Równoległa aktualizacja pozycji i kolizji ze ścianami ===
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < w->count; i++)
    {
        w->p[i].pos.x += w->p[i].vel.x * dt;
        w->p[i].pos.y += w->p[i].vel.y * dt;

        // Wall collision - lepsze: clamp pozycję + odbij prędkość
        // Lewa/prawa ściana
        if (w->p[i].pos.x - w->p[i].radius < 0)
        {
            w->p[i].pos.x = w->p[i].radius;  // Clamp do granic
            w->p[i].vel.x = -w->p[i].vel.x;  // Odbij prędkość
        }
        else if (w->p[i].pos.x + w->p[i].radius > w->width)
        {
            w->p[i].pos.x = w->width - w->p[i].radius;
            w->p[i].vel.x = -w->p[i].vel.x;
        }

        // Góra/dół ściana
        if (w->p[i].pos.y - w->p[i].radius < 0)
        {
            w->p[i].pos.y = w->p[i].radius;
            w->p[i].vel.y = -w->p[i].vel.y;
        }
        else if (w->p[i].pos.y + w->p[i].radius > w->height)
        {
            w->p[i].pos.y = w->height - w->p[i].radius;
            w->p[i].vel.y = -w->p[i].vel.y;
        }
    }

    // === OPENMP: Równoległa obsługa kolizji między cząstkami ===
    // Transformujemy zagnieżdżone pętle w jedną, aby OpenMP mogła zrównoleglić
    // Para (i,j) gdzie j > i → indeks liniowy: idx = i*count - i*(i+1)/2 + (j-i-1)
    int total_pairs = w->count * (w->count - 1) / 2;
    
    #pragma omp parallel for schedule(dynamic, 16)
    for (int pair_idx = 0; pair_idx < total_pairs; pair_idx++)
    {
        // Konwertuj indeks liniowy na (i, j)
        int i = 0;
        int temp = pair_idx;
        while (temp >= w->count - i - 1)
        {
            temp -= (w->count - i - 1);
            i++;
        }
        int j = i + 1 + temp;

        float dx = w->p[j].pos.x - w->p[i].pos.x;
        float dy = w->p[j].pos.y - w->p[i].pos.y;
        float dist2 = dx * dx + dy * dy;
        float radius_sum = w->p[i].radius + w->p[j].radius;

        if (dist2 < radius_sum * radius_sum)
        {
            // Simple elastic collision response
            float dist = sqrtf(dist2);
            float overlap = radius_sum - dist;

            // Move particles apart
            float nx = dx / dist;
            float ny = dy / dist;
            
            // === CRITICAL SECTION: Atomowe operacje na pozycjach cząstek ===
            #pragma omp critical
            {
                w->p[i].pos.x -= nx * overlap * 0.5f;
                w->p[i].pos.y -= ny * overlap * 0.5f;
                w->p[j].pos.x += nx * overlap * 0.5f;
                w->p[j].pos.y += ny * overlap * 0.5f;

                // Reflect velocities
                float relative_velocity_x = w->p[j].vel.x - w->p[i].vel.x;
                float relative_velocity_y = w->p[j].vel.y - w->p[i].vel.y;
                float velocity_along_normal = relative_velocity_x * nx + relative_velocity_y * ny;

                if (velocity_along_normal <= 0)  // Changed condition logic
                {
                    // Calculate restitution
                    float bounce_factor = 1.0f; // Perfectly elastic
                    float impulse = -(1 + bounce_factor) * velocity_along_normal;
                    impulse /= (1 / w->p[i].mass + 1 / w->p[j].mass);

                    w->p[i].vel.x -= impulse / w->p[i].mass * nx;
                    w->p[i].vel.y -= impulse / w->p[i].mass * ny;
                    w->p[j].vel.x += impulse / w->p[j].mass * nx;
                    w->p[j].vel.y += impulse / w->p[j].mass * ny;
                }
            }
        }
    }
}