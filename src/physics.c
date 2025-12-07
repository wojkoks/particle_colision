#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "physics.h"

#ifdef _OPENMP
#include <omp.h>
#endif

static inline void resolve_stick_circle(Particle *owner, Particle *target, float eps_dist)
{
    if (!owner->has_stick || owner->stick_len <= 0.0f)
        return;

    float nx = cosf(owner->stick_angle);
    float ny = sinf(owner->stick_angle);
    float sx0 = owner->pos.x + nx * owner->radius;
    float sy0 = owner->pos.y + ny * owner->radius;
    float sx1 = sx0 + nx * owner->stick_len;
    float sy1 = sy0 + ny * owner->stick_len;

    float vx = sx1 - sx0;
    float vy = sy1 - sy0;
    float seg_len2 = vx * vx + vy * vy;
    if (seg_len2 < eps_dist) return;
    float tx = target->pos.x - sx0;
    float ty = target->pos.y - sy0;
    float t = (tx * vx + ty * vy) / seg_len2;
    if (t < 0.0f) t = 0.0f;
    if (t > 1.0f) t = 1.0f;
    float closest_x = sx0 + t * vx;
    float closest_y = sy0 + t * vy;

    float dx = target->pos.x - closest_x;
    float dy = target->pos.y - closest_y;
    float dist2 = dx * dx + dy * dy;
    float rad = target->radius;
    if (dist2 >= rad * rad) return;
    float dist = sqrtf(dist2 + eps_dist);
    float overlap = rad - dist;
    float cnx = dx / dist;
    float cny = dy / dist;

    target->pos.x += cnx * overlap;
    target->pos.y += cny * overlap;
    owner->pos.x -= cnx * overlap * 0.25f;
    owner->pos.y -= cny * overlap * 0.25f;

    float rvx = target->vel.x - owner->vel.x;
    float rvy = target->vel.y - owner->vel.y;
    float vel_along_normal = rvx * cnx + rvy * cny;
    if (vel_along_normal > 0) return;
    float restitution = 1.0f;
    float inv_mass_owner = 1.0f / owner->mass;
    float inv_mass_target = 1.0f / target->mass;
    float impulse = -(1.0f + restitution) * vel_along_normal;
    impulse /= (inv_mass_owner + inv_mass_target);
    owner->vel.x -= impulse * inv_mass_owner * cnx;
    owner->vel.y -= impulse * inv_mass_owner * cny;
    target->vel.x += impulse * inv_mass_target * cnx;
    target->vel.y += impulse * inv_mass_target * cny;
}

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
        w->p[i].has_stick = (i % 10 == 0) ? 1 : 0; // co 10-ta cząstka ma kij
        w->p[i].stick_len = w->p[i].has_stick ? w->p[i].radius * 2.0f : 0.0f;
        w->p[i].stick_angle = 0.0f;
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
    printf("[INIT] Velocity stats: avg=%.2f min=%.2f max=%.2f (expected_max=%.2f)\n",
           vel_sum / count, vel_min, vel_max, max_speed * 1.414f);
}

void world_free(World *w)
{
    free(w->p);
    w->p = NULL;

#ifdef USE_CUDA
    world_cuda_release();
#endif
}

void world_reset(World *w, unsigned seed)
{
    int width = w->width;
    int height = w->height;
    int count = w->count;
    float max_speed = 120.0f;
    int rmin = 3;
    int rmax = 8;

    world_free(w);

    world_init(w, width, height, count, seed, rmin, rmax, max_speed);
}

#ifndef USE_CUDA
void world_step(World *w, float dt)
{
    const float eps_vel = 1e-4f;
    const float eps_dist = 1e-5f;

    // === OPENMP: Równoległa aktualizacja pozycji i kolizji ze ścianami ===
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < w->count; i++)
    {
        w->p[i].pos.x += w->p[i].vel.x * dt;
        w->p[i].pos.y += w->p[i].vel.y * dt;

        // aktualizacja kąta kija: kierunek prędkości
        float speed2 = w->p[i].vel.x * w->p[i].vel.x + w->p[i].vel.y * w->p[i].vel.y;
        if (w->p[i].has_stick && speed2 > eps_vel)
        {
            w->p[i].stick_angle = atan2f(w->p[i].vel.y, w->p[i].vel.x);
        }

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

        #pragma omp critical
        {
            if (dist2 < radius_sum * radius_sum)
            {
                // Simple elastic collision response
                float dist = sqrtf(dist2);
                float overlap = radius_sum - dist;

                // Move particles apart
                float nx = dx / dist;
                float ny = dy / dist;

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

            // kij vs kula w obu kierunkach
            resolve_stick_circle(&w->p[i], &w->p[j], eps_dist);
            resolve_stick_circle(&w->p[j], &w->p[i], eps_dist);
        }
    }
}
#endif
