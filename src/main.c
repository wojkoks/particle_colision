#include <SDL2/SDL.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include "physics.h"
#include "draw.h"
#include <omp.h>

int main()
{
    int width = 800;
    int height = 600;
    int particle_count = 120;
    float max_speed = 120.0f;
    int radius_min = 3;
    int radius_max = 8;
    unsigned seed = (unsigned)time(NULL);

<<<<<<< HEAD
    SDL_Init(SDL_INIT_VIDEO);
=======
    printf("OpenMP threads: %d\n", omp_get_max_threads());

    if (SDL_Init(SDL_INIT_VIDEO | SDL_INIT_TIMER) != 0)
    {
        fprintf(stderr, "SDL_Init Error: %s\n", SDL_GetError());
        return 1;
    }

>>>>>>> upstream/main
    SDL_Window *win = SDL_CreateWindow("Symulacja cząstek 2D", SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED, width, height, SDL_WINDOW_SHOWN);
    SDL_Renderer *ren = SDL_CreateRenderer(win, -1, SDL_RENDERER_ACCELERATED);

    World world;
    world_init(&world, width, height, particle_count, seed, radius_min, radius_max, max_speed);
    int paused = 0;

    int quit = 0;
    SDL_Event event;
    while (!quit)
    {
        while (SDL_PollEvent(&event))
        {
            if (event.type == SDL_QUIT)
            {
                quit = 1;
            }
            if (event.type != SDL_KEYDOWN)
            {
                continue;
            }

            switch (event.key.keysym.sym)
            {
            case SDLK_ESCAPE:
                quit = 1;
                break;
            case SDLK_SPACE:
                paused = !paused;
                break;
            case SDLK_r:
                world_reset(&world, (unsigned)time(NULL));
                break;
            case SDLK_UP:
                if (world.count < 5000)
                {
                    world_free(&world);
                    world_init(&world, width, height, world.count + 10,
                               (unsigned)time(NULL), radius_min, radius_max, max_speed);
                }
                break;
            case SDLK_DOWN:
                if (world.count > 20)
                {
                    world_free(&world);
                    world_init(&world, width, height, world.count - 10,
                               (unsigned)time(NULL), radius_min, radius_max, max_speed);
                }
                break;
            }
        }

        if (!paused)
        {
            world_step(&world, 1.0f / 120.0f);
        }

        SDL_SetRenderDrawColor(ren, 0, 0, 0, 255);
        SDL_RenderClear(ren);
        for (int i = 0; i < world.count; i++)
        {
            Particle *p = &world.p[i];
            SDL_SetRenderDrawColor(ren, p->r, p->g, p->b, 255);
            draw_filled_circle(ren, (int)p->pos.x, (int)p->pos.y, (int)p->radius);
        }
        SDL_RenderPresent(ren);
    }

    // Clean up
    world_free(&world);
    SDL_DestroyRenderer(ren);
    SDL_DestroyWindow(win);
    SDL_Quit();
    return 0;
}
