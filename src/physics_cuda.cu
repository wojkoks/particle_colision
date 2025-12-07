#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cuda_runtime.h>
#include "physics.h"

#define CUDA_CHECK(call)                                                       \
    do {                                                                       \
        cudaError_t err = (call);                                              \
        if (err != cudaSuccess) {                                              \
            fprintf(stderr, "CUDA error %s:%d: %s\n", __FILE__, __LINE__,      \
                    cudaGetErrorString(err));                                  \
            exit(EXIT_FAILURE);                                                \
        }                                                                      \
    } while (0)

static Particle *d_particles = NULL;
static int d_capacity = 0;
static int cuda_ready = 0;

__global__ void integrate_positions(Particle *p, int count, int width, int height, float dt)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= count) return;

    Particle *pi = &p[idx];

    float x = pi->pos.x + pi->vel.x * dt;
    float y = pi->pos.y + pi->vel.y * dt;

    if (x - pi->radius < 0.0f)
    {
        x = pi->radius;
        pi->vel.x = -pi->vel.x;
    }
    else if (x + pi->radius > (float)width)
    {
        x = (float)width - pi->radius;
        pi->vel.x = -pi->vel.x;
    }

    if (y - pi->radius < 0.0f)
    {
        y = pi->radius;
        pi->vel.y = -pi->vel.y;
    }
    else if (y + pi->radius > (float)height)
    {
        y = (float)height - pi->radius;
        pi->vel.y = -pi->vel.y;
    }

    pi->pos.x = x;
    pi->pos.y = y;
}

__global__ void collide_pairs(Particle *p, int count)
{
    int pair_idx = blockIdx.x * blockDim.x + threadIdx.x;
    int total_pairs = count * (count - 1) / 2;
    if (pair_idx >= total_pairs) return;

    int i = 0;
    int temp = pair_idx;
    while (temp >= count - i - 1)
    {
        temp -= (count - i - 1);
        i++;
    }
    int j = i + 1 + temp;

    float dx = p[j].pos.x - p[i].pos.x;
    float dy = p[j].pos.y - p[i].pos.y;
    float dist2 = dx * dx + dy * dy;
    float radius_sum = p[i].radius + p[j].radius;

    if (dist2 < radius_sum * radius_sum && dist2 > 1e-6f)
    {
        float dist = sqrtf(dist2);
        float overlap = radius_sum - dist;

        float nx = dx / dist;
        float ny = dy / dist;

        atomicAdd(&p[i].pos.x, -nx * overlap * 0.5f);
        atomicAdd(&p[i].pos.y, -ny * overlap * 0.5f);
        atomicAdd(&p[j].pos.x, nx * overlap * 0.5f);
        atomicAdd(&p[j].pos.y, ny * overlap * 0.5f);

        float relative_velocity_x = p[j].vel.x - p[i].vel.x;
        float relative_velocity_y = p[j].vel.y - p[i].vel.y;
        float velocity_along_normal = relative_velocity_x * nx + relative_velocity_y * ny;

        if (velocity_along_normal <= 0.0f)
        {
            float bounce_factor = 1.0f;
            float impulse = -(1.0f + bounce_factor) * velocity_along_normal;
            impulse /= (1.0f / p[i].mass + 1.0f / p[j].mass);

            float impulse_i = impulse / p[i].mass;
            float impulse_j = impulse / p[j].mass;

            atomicAdd(&p[i].vel.x, -impulse_i * nx);
            atomicAdd(&p[i].vel.y, -impulse_i * ny);
            atomicAdd(&p[j].vel.x, impulse_j * nx);
            atomicAdd(&p[j].vel.y, impulse_j * ny);
        }
    }
}

static void ensure_device(void)
{
    if (cuda_ready) return;

    int device_count = 0;
    CUDA_CHECK(cudaGetDeviceCount(&device_count));
    if (device_count == 0)
    {
        fprintf(stderr, "CUDA: no compatible devices found\n");
        exit(EXIT_FAILURE);
    }

    CUDA_CHECK(cudaSetDevice(0));
    cuda_ready = 1;
}

static void ensure_capacity(int count)
{
    if (count <= d_capacity) return;

    if (d_particles != NULL)
    {
        CUDA_CHECK(cudaFree(d_particles));
    }

    CUDA_CHECK(cudaMalloc(&d_particles, count * sizeof(Particle)));
    d_capacity = count;
}

extern "C" void world_cuda_release(void)
{
    if (d_particles != NULL)
    {
        CUDA_CHECK(cudaFree(d_particles));
    }
    d_particles = NULL;
    d_capacity = 0;
    cuda_ready = 0;
}

extern "C" void world_step(World *w, float dt)
{
    if (w->count <= 0) return;

    ensure_device();
    ensure_capacity(w->count);

    size_t bytes = (size_t)w->count * sizeof(Particle);
    CUDA_CHECK(cudaMemcpy(d_particles, w->p, bytes, cudaMemcpyHostToDevice));

    int threads = 256;
    int grid_particles = (w->count + threads - 1) / threads;
    integrate_positions<<<grid_particles, threads>>>(d_particles, w->count, w->width, w->height, dt);

    int total_pairs = w->count * (w->count - 1) / 2;
    if (total_pairs > 0)
    {
        int grid_pairs = (total_pairs + threads - 1) / threads;
        collide_pairs<<<grid_pairs, threads>>>(d_particles, w->count);
    }

    CUDA_CHECK(cudaPeekAtLastError());
    CUDA_CHECK(cudaDeviceSynchronize());
    CUDA_CHECK(cudaMemcpy(w->p, d_particles, bytes, cudaMemcpyDeviceToHost));
}
