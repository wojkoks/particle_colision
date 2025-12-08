#pragma once
#include <stdint.h>

#ifdef __cplusplus
extern "C"
{
#endif

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
        uint8_t r, g, b;
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

#ifdef USE_CUDA
    void world_cuda_release(void);
#endif

#ifdef __cplusplus
}
#endif
