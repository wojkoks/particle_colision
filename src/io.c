#include <stdio.h>
#include "io.h"

void save_frame_ppm(ParticleArray *pa, int frame_number) {
    char filename[64];
    snprintf(filename, sizeof(filename), "data/frame_%04d.ppm", frame_number);

    FILE *f = fopen(filename, "w");
    if (!f) return;

    fprintf(f, "P3\n%d %d\n255\n", WIDTH, HEIGHT);

    for (int y = 0; y < HEIGHT; y++) {
        for (int x = 0; x < WIDTH; x++) {
            int color = 0;
            for (int i = 0; i < pa->count; i++) {
                int px = (int)pa->array[i].x;
                int py = (int)pa->array[i].y;
                if (px == x && py == y) {
                    color = pa->array[i].color;
                    break;
                }
            }
            fprintf(f, "%d %d %d ", color, 0, 255 - color);
        }
        fprintf(f, "\n");
    }
    fclose(f);
}

void print_fps(double fps) {
    printf("FPS: %.2f\n", fps);
}