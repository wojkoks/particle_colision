#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "particles.h"
#include "physics.h"
#include "io.h"

int main(int argc, char *argv[]) {
    int N = 10;
    ParticleArray pa;
    init_particle_array(&pa, N);

    const int steps = 1000;
    clock_t start = clock();

    for (int step = 0; step < steps; step++) {
        update_positions(&pa);
        handle_collisions(&pa);
        if (step % 10 == 0)
            save_frame_ppm(&pa, step);
    }

    double time_spent = (double)(clock() - start) / CLOCKS_PER_SEC;
    print_fps(steps / time_spent);

    free_particle_array(&pa);
    return 0;
}