#pragma once
#include <SDL2/SDL.h>
#include <mpi.h>
#include "physics.h"

// === MPI Domain Decomposition - 2D Cartesian Grid ===
// Plansza podzielona na siatkę 2x2 (lub większą)
// Każdy proces zarządza blokiem i wymienia cząstkami z sąsiadami

typedef struct
{
    int rank;           // ID procesu (0, 1, 2, ...)
    int world_size;     // Liczba procesów razem
    
    // === Topologia 2D (Cartesian) ===
    MPI_Comm cart_comm;   // Communicator z topologią kartezjańską
    int dims[2];          // Wymiary siatki [rows, cols]
    int coords[2];        // Moje współrzędne [row, col]
    int nbr_up, nbr_down, nbr_left, nbr_right;  // Sąsiedzi główni
    int nbr_ul, nbr_ur, nbr_dl, nbr_dr;         // Sąsiedzi diagonalni
    
    // === Wymiary ===
    int global_width;     // Całkowita szerokość (800)
    int global_height;    // Całkowita wysokość (600)
    
    // === Mój blok (domena tego procesu) ===
    int block_x0, block_y0;       // Górny-lewy róg
    int block_width, block_height; // Rozmiar
    
    // === Cząstki ===
    Particle *particles;       // Tablica cząstek w moim bloku
    int particle_count;        // Liczba cząstek
    int particle_capacity;     // Alokowana pojemność
    
    // === Ghost cells - kopie particles z sąsiednich bloków ===
    Particle *ghost_particles;  // Kopie particles ze sąsiadów dla kolizji
    int ghost_count;            // Liczba ghost particles
    
    // === Bufory wymiany (do każdego sąsiada) ===
    Particle *send_buf[8];     // Bufory wysłania [up, down, left, right, ul, ur, dl, dr]
    int send_count[8];
    Particle *recv_buf[8];     // Bufory odbierania
    int recv_count[8];
    int buf_capacity;
    
} WorldMPI;

// === Funkcje ===
void world_mpi_init(WorldMPI *w, int global_width, int global_height,
                    int total_particles, unsigned seed, int rmin, int rmax, float max_speed);
void world_mpi_free(WorldMPI *w);
void world_mpi_step(WorldMPI *w, float dt);
void world_mpi_gather_to_root(WorldMPI *w, Particle **all_particles, int *total_count);
