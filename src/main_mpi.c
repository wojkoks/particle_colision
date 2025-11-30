#include <SDL2/SDL.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#include <sys/types.h>
#include <mpi.h>
#include "physics_mpi.h"
#include "draw.h"

int main(int argc, char *argv[])
{
    // === MPI Initialization ===
    MPI_Init(&argc, &argv);

    int rank, world_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // === Parametry symulacji ===
    int width = 800;
    int height = 600;
    int particle_count = 120;
    float max_speed = 120.0f;
    int radius_min = 3;
    int radius_max = 8;
    unsigned seed = (unsigned)time(NULL);

    // === Inicjalizuj świat MPI (wszyscy procesy) ===
    WorldMPI world;
    world_mpi_init(&world, width, height, particle_count, seed, radius_min, radius_max, max_speed);

    // ========================================================================
    // === PROCESY 1+ : Pracują w tle bez SDL ===
    // ========================================================================
    if (rank != 0)
    {
        int quit = 0;
        int paused = 0;
        int frame = 0;

        while (!quit)
        {
            // Czekaj na sygnały z procesu 0
            MPI_Bcast(&quit, 1, MPI_INT, 0, world.cart_comm);
            if (quit) break;

            MPI_Bcast(&paused, 1, MPI_INT, 0, world.cart_comm);

            // Robi krok symulacji
            if (!paused)
            {
                world_mpi_step(&world, 1.0f / 120.0f);
            }

            // Zbierz co 2 klatki (mniej transferu danych)
            if (frame % 2 == 0)
            {
                Particle *dummy = NULL;
                int dummy_count = 0;
                world_mpi_gather_to_root(&world, &dummy, &dummy_count);
                if (dummy != NULL) free(dummy);
            }

            // Debug output every 120 frames
            frame++;
            if (frame % 120 == 0)
            {
                printf("[MPI] Rank %d | Frame %d | Particles: %d\n",
                       rank, frame, world.particle_count);
            }

            // Lekki delay aby nie busy-wait
            usleep(100);
        }

        // Cleanup dla proc != 0
        world_mpi_free(&world);
        MPI_Finalize();
        return 0;
    }

    // ========================================================================
    // === PROCES 0 : SDL + zarządzanie symulacją ===
    // ========================================================================

    // Inicjalizuj SDL
    if (SDL_Init(SDL_INIT_VIDEO | SDL_INIT_TIMER) != 0)
    {
        fprintf(stderr, "SDL_Init Error: %s\n", SDL_GetError());
        MPI_Finalize();
        return 1;
    }

    SDL_Window *win = SDL_CreateWindow("Symulacja cząstek 2D (MPI)",
                                       SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED,
                                       width, height, SDL_WINDOW_SHOWN);
    if (win == NULL)
    {
        fprintf(stderr, "SDL_CreateWindow Error: %s\n", SDL_GetError());
        SDL_Quit();
        MPI_Finalize();
        return 1;
    }

    SDL_Renderer *ren = SDL_CreateRenderer(win, -1, SDL_RENDERER_ACCELERATED);
    if (ren == NULL)
    {
        fprintf(stderr, "SDL_CreateRenderer Error: %s\n", SDL_GetError());
        SDL_DestroyWindow(win);
        SDL_Quit();
        MPI_Finalize();
        return 1;
    }

    printf("[MPI] Process 0: Initialized %d processes\n", world_size);

    int paused = 0;
    int quit = 0;
    int frame = 0;

    // === Main loop - Proces 0 ===
    while (!quit)
    {
        // Obsługa eventów SDL (non-blocking)
        SDL_Event event;
        while (SDL_PollEvent(&event))
        {
            if (event.type == SDL_QUIT)
            {
                quit = 1;
            }
            if (event.type == SDL_KEYDOWN)
            {
                switch (event.key.keysym.sym)
                {
                case SDLK_ESCAPE:
                    quit = 1;
                    break;
                case SDLK_SPACE:
                    paused = !paused;
                    printf("[MPI] Paused: %s\n", paused ? "YES" : "NO");
                    break;
                }
            }
        }

        // === BROADCAST do wszystkich procesów (używamy cartesian comm) ===
        MPI_Bcast(&quit, 1, MPI_INT, 0, world.cart_comm);
        if (quit) break;

        MPI_Bcast(&paused, 1, MPI_INT, 0, world.cart_comm);

        // === Krok symulacji (wszyscy procesy) ===
        if (!paused)
        {
            world_mpi_step(&world, 1.0f / 120.0f);
        }

        // === BRAK MPI_Barrier - procesy nie czekają na siebie ===
        // (MPI_Barrier usuwaliliśmy - zbyt dużo opóźnia)

        // === Zbierz cząstki do rysowania (co 2 klatki) ===
        if (frame % 2 == 0)
        {
            Particle *all_particles = NULL;
            int total_particles = 0;
            world_mpi_gather_to_root(&world, &all_particles, &total_particles);

            // === Rysuj (tylko proces 0) ===
            if (rank == 0)
            {
                SDL_SetRenderDrawColor(ren, 0, 0, 0, 255);
                SDL_RenderClear(ren);

                for (int i = 0; i < total_particles; i++)
                {
                    Particle *p = &all_particles[i];
                    SDL_SetRenderDrawColor(ren, p->r, p->g, p->b, 255);
                    draw_filled_circle(ren, (int)p->pos.x, (int)p->pos.y, (int)p->radius);
                }

                // Rysuj linie gridu (podziały między procesami)
                draw_grid(ren, width, height, world.dims[1], world.dims[0]);

                SDL_RenderPresent(ren);

                // FPS control - zmniejszone z 8ms na 2ms
                SDL_Delay(2);

                // Debug info
                if (frame % 120 == 0)
                {
                    printf("[MPI] Frame %d | Total: %d | Rank 0: %d\n", 
                           frame, total_particles, world.particle_count);
                }
            }

            if (all_particles != NULL)
                free(all_particles);
        }
        else
        {
            // Na parzystych klatkach: pomiń gather, ale nie blokuj procesów
            if (rank == 0)
                SDL_Delay(1);  // Lekki delay by nie busy-loop
        }

        frame++;

    }

    // === Cleanup ===
    world_mpi_free(&world);
    SDL_DestroyRenderer(ren);
    SDL_DestroyWindow(win);
    SDL_Quit();

    MPI_Finalize();
    return 0;
}
