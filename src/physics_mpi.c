#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include <omp.h>
#include "physics_mpi.h"

// === MPI Datatype dla Particle ===
static MPI_Datatype PARTICLE_MPI_TYPE = MPI_DATATYPE_NULL;

void create_particle_mpi_type()
{
    if (PARTICLE_MPI_TYPE != MPI_DATATYPE_NULL)
        return;
    // Use contiguous bytes - simpler and safer
    MPI_Type_contiguous(sizeof(Particle), MPI_BYTE, &PARTICLE_MPI_TYPE);
    MPI_Type_commit(&PARTICLE_MPI_TYPE);
}

// === Topologia 2D (Cartesian Grid) ===
void setup_2d_cartesian(WorldMPI *w)
{
    // Oblicz wymiary siatki - najpierw spróbuj rozsądnie
    w->dims[0] = 1;
    w->dims[1] = w->world_size;
    
    // Spróbuj znaleźć bardziej kwadratową siatkę
    for (int i = 1; i * i <= w->world_size; i++)
    {
        if (w->world_size % i == 0)
        {
            w->dims[0] = i;
            w->dims[1] = w->world_size / i;
        }
    }
    
    // Upewnij się że wymiary są prawidłowe
    if (w->dims[0] * w->dims[1] != w->world_size)
    {
        w->dims[0] = 1;
        w->dims[1] = w->world_size;
    }
    
    // Utwórz topologię kartezjańską (bez periodyczności)
    int periods[2] = {0, 0};  // Nie-periodyczna (nie zawijamy na krawędziach)
    MPI_Cart_create(MPI_COMM_WORLD, 2, w->dims, periods, 1, &w->cart_comm);
    
    // Odczytaj moje współrzędne
    MPI_Cart_coords(w->cart_comm, w->rank, 2, w->coords);
    
    // Pobierz rangi sąsiadów (4-kierunkowych)
    MPI_Cart_shift(w->cart_comm, 0, 1, &w->nbr_up,    &w->nbr_down);   // Rzęd
    MPI_Cart_shift(w->cart_comm, 1, 1, &w->nbr_left,  &w->nbr_right);  // Kolumna
    
    // Diagonalni sąsiedzi (musimy obliczyć ręcznie)
    int coords_tmp[2];
    int rank_tmp;
    
    // Górny-lewy
    coords_tmp[0] = w->coords[0] - 1; coords_tmp[1] = w->coords[1] - 1;
    w->nbr_ul = (coords_tmp[0] >= 0 && coords_tmp[1] >= 0) 
        ? (MPI_Cart_rank(w->cart_comm, coords_tmp, &rank_tmp), rank_tmp) : MPI_PROC_NULL;
    
    // Górny-prawy
    coords_tmp[0] = w->coords[0] - 1; coords_tmp[1] = w->coords[1] + 1;
    w->nbr_ur = (coords_tmp[0] >= 0 && coords_tmp[1] < w->dims[1]) 
        ? (MPI_Cart_rank(w->cart_comm, coords_tmp, &rank_tmp), rank_tmp) : MPI_PROC_NULL;
    
    // Dolny-lewy
    coords_tmp[0] = w->coords[0] + 1; coords_tmp[1] = w->coords[1] - 1;
    w->nbr_dl = (coords_tmp[0] < w->dims[0] && coords_tmp[1] >= 0) 
        ? (MPI_Cart_rank(w->cart_comm, coords_tmp, &rank_tmp), rank_tmp) : MPI_PROC_NULL;
    
    // Dolny-prawy
    coords_tmp[0] = w->coords[0] + 1; coords_tmp[1] = w->coords[1] + 1;
    w->nbr_dr = (coords_tmp[0] < w->dims[0] && coords_tmp[1] < w->dims[1]) 
        ? (MPI_Cart_rank(w->cart_comm, coords_tmp, &rank_tmp), rank_tmp) : MPI_PROC_NULL;
    
    // Oblicz wymiary każdego bloku
    int block_width  = w->global_width / w->dims[1];   // Kolumny
    int block_height = w->global_height / w->dims[0];  // Rzędy
    
    w->block_x0 = w->coords[1] * block_width;
    w->block_y0 = w->coords[0] * block_height;
    w->block_width  = (w->coords[1] == w->dims[1] - 1) 
        ? (w->global_width  - w->block_x0) : block_width;
    w->block_height = (w->coords[0] == w->dims[0] - 1) 
        ? (w->global_height - w->block_y0) : block_height;
}

void world_mpi_init(WorldMPI *w, int global_width, int global_height,
                    int total_particles, unsigned seed, int rmin, int rmax, float max_speed)
{
    MPI_Comm_rank(MPI_COMM_WORLD, &w->rank);
    MPI_Comm_size(MPI_COMM_WORLD, &w->world_size);

    create_particle_mpi_type();

    w->global_width = global_width;
    w->global_height = global_height;

    // Setup topologia 2D
    setup_2d_cartesian(w);

    // === Alokacja cząstek ===
    // Każdy proces dostaje ~równo cząstek
    int particles_per_process = total_particles / w->world_size;
    int extra = total_particles % w->world_size;

    w->particle_capacity = particles_per_process + (w->rank == 0 ? extra : 0) + 50;
    w->particles = (Particle *)malloc(w->particle_capacity * sizeof(Particle));
    w->particle_count = (w->rank == 0) ? (particles_per_process + extra) : particles_per_process;

    // === Inicjalizacja cząstek na procesie 0, rozsyłanie na inne ===
    if (w->rank == 0)
    {
        srand(seed);
        // Rank 0 tworzy particles dla WSZYSTKICH procesów
        // ale każdy process dostaje particles tylko w jego obszarze
        
        Particle temp_all_particles[total_particles];
        int particle_idx = 0;
        
        // Dla każdego procesu, utwórz particles w jego bloku
        for (int proc_rank = 0; proc_rank < w->world_size; proc_rank++)
        {
            int proc_particles = (proc_rank == w->world_size - 1) 
                ? (total_particles - particle_idx) 
                : (total_particles / w->world_size);
            
            // Oblicz blok tego procesu (uproszczone - zakładając 2D grid)
            int block_col = proc_rank % w->dims[1];
            int block_row = proc_rank / w->dims[1];
            
            int block_x0 = (global_width / w->dims[1]) * block_col;
            int block_y0 = (global_height / w->dims[0]) * block_row;
            int block_width = (block_col == w->dims[1] - 1) 
                ? (global_width - block_x0) 
                : (global_width / w->dims[1]);
            int block_height = (block_row == w->dims[0] - 1) 
                ? (global_height - block_y0) 
                : (global_height / w->dims[0]);
            
            // Utwórz particles w tym bloku
            for (int i = 0; i < proc_particles; i++)
            {
                if (particle_idx >= total_particles) break;
                
                temp_all_particles[particle_idx].radius = rand() % (rmax - rmin + 1) + rmin;
                
                // Pozycja w obrębie bloku tego procesu!
                float margin = (float)temp_all_particles[particle_idx].radius;
                temp_all_particles[particle_idx].pos.x = block_x0 + margin + 
                    (float)(rand()) / RAND_MAX * (block_width - 2 * margin);
                temp_all_particles[particle_idx].pos.y = block_y0 + margin + 
                    (float)(rand()) / RAND_MAX * (block_height - 2 * margin);
                
                temp_all_particles[particle_idx].vel.x = ((float)rand() / RAND_MAX) * max_speed * (rand() % 2 == 0 ? 1 : -1);
                temp_all_particles[particle_idx].vel.y = ((float)rand() / RAND_MAX) * max_speed * (rand() % 2 == 0 ? 1 : -1);
                temp_all_particles[particle_idx].mass = 3.14f * temp_all_particles[particle_idx].radius * temp_all_particles[particle_idx].radius;
                temp_all_particles[particle_idx].r = rand() % 256;
                temp_all_particles[particle_idx].g = rand() % 256;
                temp_all_particles[particle_idx].b = rand() % 256;
                
                particle_idx++;
            }
        }
        
        // Rozpakuj do swoich particles
        int my_count = total_particles / w->world_size;
        for (int i = 0; i < my_count; i++)
        {
            w->particles[i] = temp_all_particles[i];
        }
        w->particle_count = my_count;
        
        // Debug: calculate velocity stats (only rank 0)
        if (w->rank == 0)
        {
            float vel_sum = 0, vel_max = 0, vel_min = 999;
            for (int i = 0; i < total_particles; i++)
            {
                float speed = sqrtf(temp_all_particles[i].vel.x * temp_all_particles[i].vel.x + 
                                   temp_all_particles[i].vel.y * temp_all_particles[i].vel.y);
                vel_sum += speed;
                if (speed > vel_max) vel_max = speed;
                if (speed < vel_min) vel_min = speed;
            }
            printf("[MPI] Velocity stats: avg=%.2f min=%.2f max=%.2f (expected_max=%.2f)\n", 
                   vel_sum / total_particles, vel_min, vel_max, max_speed * 1.414f);
        }

        int idx = my_count;
        for (int rank = 1; rank < w->world_size; rank++)
        {
            int count = (rank == w->world_size - 1) 
                ? (total_particles - idx) 
                : (total_particles / w->world_size);
            MPI_Send(&temp_all_particles[idx], count, PARTICLE_MPI_TYPE, rank, 0, MPI_COMM_WORLD);
            idx += count;
        }
    }
    else
    {
        MPI_Recv(w->particles, w->particle_count, PARTICLE_MPI_TYPE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    // === Bufory wymiany (dla 8 sąsiadów) ===
    // Make buffer larger to avoid losing particles
    w->buf_capacity = w->particle_capacity;  // Full capacity per buffer
    for (int i = 0; i < 8; i++)
    {
        w->send_buf[i] = (Particle *)malloc(w->buf_capacity * sizeof(Particle));
        w->recv_buf[i] = (Particle *)malloc(w->buf_capacity * sizeof(Particle));
        w->send_count[i] = 0;
        w->recv_count[i] = 0;
    }
    
    // Ghost cells - kopie particles z sąsiednich bloków do kolizji
    w->ghost_particles = (Particle *)malloc(w->particle_capacity * 4 * sizeof(Particle));
    w->ghost_count = 0;

    if (w->rank == 0)
    {
        printf("[MPI] Grid: %dx%d, Total processes: %d\n", w->dims[0], w->dims[1], w->world_size);
        printf("[MPI] Rank 0: block=[%d,%d] size=%dx%d, received %d particles\n",
               w->block_x0, w->block_y0, w->block_width, w->block_height, w->particle_count);
    }
    else
    {
        printf("[MPI] Rank %d: block=[%d,%d] size=%dx%d, received %d particles\n",
               w->rank, w->block_x0, w->block_y0, w->block_width, w->block_height, w->particle_count);
    }
}

void world_mpi_free(WorldMPI *w)
{
    free(w->particles);
    free(w->ghost_particles);
    for (int i = 0; i < 8; i++)
    {
        free(w->send_buf[i]);
        free(w->recv_buf[i]);
    }
    w->particles = NULL;
    w->ghost_particles = NULL;

    if (PARTICLE_MPI_TYPE != MPI_DATATYPE_NULL)
    {
        MPI_Type_free(&PARTICLE_MPI_TYPE);
        PARTICLE_MPI_TYPE = MPI_DATATYPE_NULL;
    }
}

void world_mpi_step(WorldMPI *w, float dt)
{
    // === Krok 1: Aktualizacja pozycji (bez OpenMP - szybciej) ===
    for (int i = 0; i < w->particle_count; i++)
    {
        w->particles[i].pos.x += w->particles[i].vel.x * dt;
        w->particles[i].pos.y += w->particles[i].vel.y * dt;

        // Globalne granice (cały świat)
        if (w->particles[i].pos.x - w->particles[i].radius < 0)
        {
            w->particles[i].pos.x = w->particles[i].radius;
            w->particles[i].vel.x = -w->particles[i].vel.x;
        }
        else if (w->particles[i].pos.x + w->particles[i].radius > w->global_width)
        {
            w->particles[i].pos.x = w->global_width - w->particles[i].radius;
            w->particles[i].vel.x = -w->particles[i].vel.x;
        }

        if (w->particles[i].pos.y - w->particles[i].radius < 0)
        {
            w->particles[i].pos.y = w->particles[i].radius;
            w->particles[i].vel.y = -w->particles[i].vel.y;
        }
        else if (w->particles[i].pos.y + w->particles[i].radius > w->global_height)
        {
            w->particles[i].pos.y = w->global_height - w->particles[i].radius;
            w->particles[i].vel.y = -w->particles[i].vel.y;
        }
    }

    // === Krok 2: Kolizje między cząstkami (lokalne, OpenMP z reduction) ===
    int total_pairs = w->particle_count * (w->particle_count - 1) / 2;

    #pragma omp parallel for schedule(dynamic)
    for (int pair_idx = 0; pair_idx < total_pairs; pair_idx++)
    {
        int i = 0;
        int temp = pair_idx;
        while (temp >= w->particle_count - i - 1)
        {
            temp -= (w->particle_count - i - 1);
            i++;
        }
        int j = i + 1 + temp;

        float dx = w->particles[j].pos.x - w->particles[i].pos.x;
        float dy = w->particles[j].pos.y - w->particles[i].pos.y;
        float dist2 = dx * dx + dy * dy;
        float radius_sum = w->particles[i].radius + w->particles[j].radius;

        if (dist2 < radius_sum * radius_sum)
        {
            float dist = sqrtf(dist2);
            float overlap = radius_sum - dist;

            float nx = dx / dist;
            float ny = dy / dist;

            // Atomowe operacje zamiast critical
            #pragma omp atomic
            w->particles[i].pos.x -= nx * overlap * 0.5f;
            #pragma omp atomic
            w->particles[i].pos.y -= ny * overlap * 0.5f;
            #pragma omp atomic
            w->particles[j].pos.x += nx * overlap * 0.5f;
            #pragma omp atomic
            w->particles[j].pos.y += ny * overlap * 0.5f;

            float relative_velocity_x = w->particles[j].vel.x - w->particles[i].vel.x;
            float relative_velocity_y = w->particles[j].vel.y - w->particles[i].vel.y;
            float velocity_along_normal = relative_velocity_x * nx + relative_velocity_y * ny;

            if (velocity_along_normal <= 0)
            {
                float bounce_factor = 1.0f;
                float impulse = -(1 + bounce_factor) * velocity_along_normal;
                impulse /= (1 / w->particles[i].mass + 1 / w->particles[j].mass);

                #pragma omp atomic
                w->particles[i].vel.x -= impulse / w->particles[i].mass * nx;
                #pragma omp atomic
                w->particles[i].vel.y -= impulse / w->particles[i].mass * ny;
                #pragma omp atomic
                w->particles[j].vel.x += impulse / w->particles[j].mass * nx;
                #pragma omp atomic
                w->particles[j].vel.y += impulse / w->particles[j].mass * ny;
            }
        }
    }

    // === Krok 2.5: Wymiana ghost cells z wszystkimi 4 sąsiadami (bez OpenMP) ===
    w->ghost_count = 0;
    MPI_Status status;
    
    Particle temp_send[w->particle_capacity];
    Particle temp_recv[w->particle_capacity];
    int temp_count = 0;
    int recv_count = 0;
    
    float ghost_margin = 100.0f;  // Margines dla ghost particles (większy dla pewności)
    
    // === Wyślij do UP (nbr_up), otrzymaj z DOŁU (nbr_down) ===
    if (w->nbr_up != MPI_PROC_NULL || w->nbr_down != MPI_PROC_NULL)
    {
        temp_count = 0;
        // Zbierz particles z górnej części bloku
        for (int i = 0; i < w->particle_count; i++)
        {
            if (w->particles[i].pos.y <= w->block_y0 + ghost_margin)
            {
                if (temp_count < w->particle_capacity)
                    temp_send[temp_count++] = w->particles[i];
            }
        }
        
        // Wymień
        recv_count = 0;
        MPI_Sendrecv(temp_send, temp_count, PARTICLE_MPI_TYPE,
                     w->nbr_up, 10,
                     temp_recv, w->particle_capacity, PARTICLE_MPI_TYPE,
                     w->nbr_down, 10,
                     w->cart_comm, &status);
        MPI_Get_count(&status, PARTICLE_MPI_TYPE, &recv_count);
        for (int i = 0; i < recv_count && w->ghost_count < w->particle_capacity; i++)
        {
            w->ghost_particles[w->ghost_count++] = temp_recv[i];
        }
    }
    
    // === Wyślij do DOWN (nbr_down), otrzymaj z GÓRY (nbr_up) ===
    if (w->nbr_down != MPI_PROC_NULL || w->nbr_up != MPI_PROC_NULL)
    {
        temp_count = 0;
        // Zbierz particles z dolnej części bloku
        for (int i = 0; i < w->particle_count; i++)
        {
            if (w->particles[i].pos.y >= w->block_y0 + w->block_height - ghost_margin)
            {
                if (temp_count < w->particle_capacity)
                    temp_send[temp_count++] = w->particles[i];
            }
        }
        
        // Wymień
        recv_count = 0;
        MPI_Sendrecv(temp_send, temp_count, PARTICLE_MPI_TYPE,
                     w->nbr_down, 11,
                     temp_recv, w->particle_capacity, PARTICLE_MPI_TYPE,
                     w->nbr_up, 11,
                     w->cart_comm, &status);
        MPI_Get_count(&status, PARTICLE_MPI_TYPE, &recv_count);
        for (int i = 0; i < recv_count && w->ghost_count < w->particle_capacity; i++)
        {
            w->ghost_particles[w->ghost_count++] = temp_recv[i];
        }
    }
    
    // === Wyślij do LEFT (nbr_left), otrzymaj z PRAWEJ (nbr_right) ===
    if (w->nbr_left != MPI_PROC_NULL || w->nbr_right != MPI_PROC_NULL)
    {
        temp_count = 0;
        // Zbierz particles z lewej części bloku
        for (int i = 0; i < w->particle_count; i++)
        {
            if (w->particles[i].pos.x <= w->block_x0 + ghost_margin)
            {
                if (temp_count < w->particle_capacity)
                    temp_send[temp_count++] = w->particles[i];
            }
        }
        
        // Wymień
        recv_count = 0;
        MPI_Sendrecv(temp_send, temp_count, PARTICLE_MPI_TYPE,
                     w->nbr_left, 12,
                     temp_recv, w->particle_capacity, PARTICLE_MPI_TYPE,
                     w->nbr_right, 12,
                     w->cart_comm, &status);
        MPI_Get_count(&status, PARTICLE_MPI_TYPE, &recv_count);
        for (int i = 0; i < recv_count && w->ghost_count < w->particle_capacity; i++)
        {
            w->ghost_particles[w->ghost_count++] = temp_recv[i];
        }
    }
    
    // === Wyślij do RIGHT (nbr_right), otrzymaj z LEWEJ (nbr_left) ===
    if (w->nbr_right != MPI_PROC_NULL || w->nbr_left != MPI_PROC_NULL)
    {
        temp_count = 0;
        // Zbierz particles z prawej części bloku
        for (int i = 0; i < w->particle_count; i++)
        {
            if (w->particles[i].pos.x >= w->block_x0 + w->block_width - ghost_margin)
            {
                if (temp_count < w->particle_capacity)
                    temp_send[temp_count++] = w->particles[i];
            }
        }
        
        // Wymień
        recv_count = 0;
        MPI_Sendrecv(temp_send, temp_count, PARTICLE_MPI_TYPE,
                     w->nbr_right, 13,
                     temp_recv, w->particle_capacity, PARTICLE_MPI_TYPE,
                     w->nbr_left, 13,
                     w->cart_comm, &status);
        MPI_Get_count(&status, PARTICLE_MPI_TYPE, &recv_count);
        for (int i = 0; i < recv_count && w->ghost_count < w->particle_capacity; i++)
        {
            w->ghost_particles[w->ghost_count++] = temp_recv[i];
        }
    }

    // === Krok 3: Kolizje z ghost particles (bez OpenMP) ===
    // Teraz sprawdzaj kolizje między własnymi i ghost particles
    for (int i = 0; i < w->particle_count; i++)
    {
        for (int j = 0; j < w->ghost_count; j++)
        {
            float dx = w->ghost_particles[j].pos.x - w->particles[i].pos.x;
            float dy = w->ghost_particles[j].pos.y - w->particles[i].pos.y;
            float dist2 = dx * dx + dy * dy;
            float radius_sum = w->particles[i].radius + w->ghost_particles[j].radius;

            if (dist2 < radius_sum * radius_sum && dist2 > 0)
            {
                float dist = sqrtf(dist2);
                float overlap = radius_sum - dist;

                float nx = dx / dist;
                float ny = dy / dist;

                // Tylko odsunięcie (nie zmieniaj ghost particles)
                w->particles[i].pos.x -= nx * overlap;
                w->particles[i].pos.y -= ny * overlap;

                float relative_velocity_x = w->ghost_particles[j].vel.x - w->particles[i].vel.x;
                float relative_velocity_y = w->ghost_particles[j].vel.y - w->particles[i].vel.y;
                float velocity_along_normal = relative_velocity_x * nx + relative_velocity_y * ny;

                if (velocity_along_normal <= 0)
                {
                    float bounce_factor = 1.0f;
                    float impulse = -(1 + bounce_factor) * velocity_along_normal;
                    impulse /= (1 / w->particles[i].mass + 1 / w->ghost_particles[j].mass);

                    w->particles[i].vel.x -= impulse / w->particles[i].mass * nx;
                    w->particles[i].vel.y -= impulse / w->particles[i].mass * ny;
                }
            }
        }
    }

    // === Krok 4: Particle migration - wyślij particles które wychodzą z bloku (bez OpenMP) ===
    // Zbierz particles do wysłania każdemu sąsiadowi
    Particle send_up[w->particle_capacity];
    Particle send_down[w->particle_capacity];
    Particle send_left[w->particle_capacity];
    Particle send_right[w->particle_capacity];
    int count_up = 0, count_down = 0, count_left = 0, count_right = 0;
    
    // Usuwamy particles które wychodzą i dodajemy do send buffers
    Particle remaining[w->particle_capacity];
    int remaining_count = 0;
    
    for (int i = 0; i < w->particle_count; i++)
    {
        Particle *p = &w->particles[i];
        int moved = 0;
        
        // Czy wychodzi na górę?
        if (p->pos.y - p->radius < w->block_y0 && w->nbr_up != MPI_PROC_NULL)
        {
            if (count_up < w->particle_capacity)
                send_up[count_up++] = *p;
            moved = 1;
        }
        // Czy wychodzi na dół?
        else if (p->pos.y + p->radius > w->block_y0 + w->block_height && w->nbr_down != MPI_PROC_NULL)
        {
            if (count_down < w->particle_capacity)
                send_down[count_down++] = *p;
            moved = 1;
        }
        // Czy wychodzi na lewo?
        else if (p->pos.x - p->radius < w->block_x0 && w->nbr_left != MPI_PROC_NULL)
        {
            if (count_left < w->particle_capacity)
                send_left[count_left++] = *p;
            moved = 1;
        }
        // Czy wychodzi na prawo?
        else if (p->pos.x + p->radius > w->block_x0 + w->block_width && w->nbr_right != MPI_PROC_NULL)
        {
            if (count_right < w->particle_capacity)
                send_right[count_right++] = *p;
            moved = 1;
        }
        
        // Jeśli nie wyszedł, zostaje w tym procesie
        if (!moved)
        {
            if (remaining_count < w->particle_capacity)
                remaining[remaining_count++] = *p;
        }
    }
    
    // Wymień particles z sąsiadami i otrzymaj ich particles
    Particle recv_up[w->particle_capacity];
    Particle recv_down[w->particle_capacity];
    Particle recv_left[w->particle_capacity];
    Particle recv_right[w->particle_capacity];
    int count_recv_up = 0, count_recv_down = 0, count_recv_left = 0, count_recv_right = 0;
    
    // UP exchange
    if (w->nbr_up != MPI_PROC_NULL || w->nbr_down != MPI_PROC_NULL)
    {
        MPI_Sendrecv(send_up, count_up, PARTICLE_MPI_TYPE,
                     w->nbr_up, 20,
                     recv_down, w->particle_capacity, PARTICLE_MPI_TYPE,
                     w->nbr_down, 20,
                     w->cart_comm, &status);
        MPI_Get_count(&status, PARTICLE_MPI_TYPE, &count_recv_down);
    }
    
    // DOWN exchange
    if (w->nbr_down != MPI_PROC_NULL || w->nbr_up != MPI_PROC_NULL)
    {
        MPI_Sendrecv(send_down, count_down, PARTICLE_MPI_TYPE,
                     w->nbr_down, 21,
                     recv_up, w->particle_capacity, PARTICLE_MPI_TYPE,
                     w->nbr_up, 21,
                     w->cart_comm, &status);
        MPI_Get_count(&status, PARTICLE_MPI_TYPE, &count_recv_up);
    }
    
    // LEFT exchange
    if (w->nbr_left != MPI_PROC_NULL || w->nbr_right != MPI_PROC_NULL)
    {
        MPI_Sendrecv(send_left, count_left, PARTICLE_MPI_TYPE,
                     w->nbr_left, 22,
                     recv_right, w->particle_capacity, PARTICLE_MPI_TYPE,
                     w->nbr_right, 22,
                     w->cart_comm, &status);
        MPI_Get_count(&status, PARTICLE_MPI_TYPE, &count_recv_right);
    }
    
    // RIGHT exchange
    if (w->nbr_right != MPI_PROC_NULL || w->nbr_left != MPI_PROC_NULL)
    {
        MPI_Sendrecv(send_right, count_right, PARTICLE_MPI_TYPE,
                     w->nbr_right, 23,
                     recv_left, w->particle_capacity, PARTICLE_MPI_TYPE,
                     w->nbr_left, 23,
                     w->cart_comm, &status);
        MPI_Get_count(&status, PARTICLE_MPI_TYPE, &count_recv_left);
    }
    
    // Aktualizuj lokalną listę particles
    w->particle_count = remaining_count;
    for (int i = 0; i < remaining_count; i++)
    {
        w->particles[i] = remaining[i];
    }
    
    // Dodaj received particles
    for (int i = 0; i < count_recv_up && w->particle_count < w->particle_capacity; i++)
    {
        w->particles[w->particle_count++] = recv_up[i];
    }
    for (int i = 0; i < count_recv_down && w->particle_count < w->particle_capacity; i++)
    {
        w->particles[w->particle_count++] = recv_down[i];
    }
    for (int i = 0; i < count_recv_left && w->particle_count < w->particle_capacity; i++)
    {
        w->particles[w->particle_count++] = recv_left[i];
    }
    for (int i = 0; i < count_recv_right && w->particle_count < w->particle_capacity; i++)
    {
        w->particles[w->particle_count++] = recv_right[i];
    }

}

void world_mpi_gather_to_root(WorldMPI *w, Particle **all_particles, int *total_count)
{
    int *counts = (int *)malloc(w->world_size * sizeof(int));
    int *displs = (int *)malloc(w->world_size * sizeof(int));

    if (w->rank == 0)
    {
        *all_particles = (Particle *)malloc(100000 * sizeof(Particle));
    }

    // === Gather - wszyscy podają liczbę swoich cząstek ===
    MPI_Gather(&w->particle_count, 1, MPI_INT, 
               counts, 1, MPI_INT, 
               0, w->cart_comm);

    // === Oblicz displacement (tylko na rank 0) ===
    if (w->rank == 0)
    {
        *total_count = 0;
        displs[0] = 0;
        for (int i = 0; i < w->world_size; i++)
        {
            *total_count += counts[i];
            if (i > 0)
                displs[i] = displs[i - 1] + counts[i - 1];
        }
    }

    // === Gatherv - zbierz wszystkie cząstki ===
    MPI_Gatherv(w->particles, w->particle_count, PARTICLE_MPI_TYPE,
                w->rank == 0 ? *all_particles : NULL,
                counts, displs, PARTICLE_MPI_TYPE,
                0, w->cart_comm);

    free(counts);
    free(displs);
}
