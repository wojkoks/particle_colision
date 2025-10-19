#include "physics.h"
#include <math.h>

void update_positions(ParticleArray *pa){
    for(int i = 0; i < pa->count; i++){
        Particle *p = &pa->array[i];
        p->x += p->vx * DT;
        p->y += p->vy * DT;
        bounce_from_walls(p);
    }
}

void bounce_from_walls(Particle *p){
    if(p->x <=0 || p->x >= WIDTH){
        p->vx = -p->vx;
    }
    if(p->y <=0 || p->y >= HEIGHT){
        p->vy = -p->vy;
    }
}

void handle_collisions(ParticleArray *pa) {
    for (int i = 0; i < pa->count; i++) {
        for (int j = i + 1; j < pa->count; j++) {
            float dx = pa->array[i].x - pa->array[j].x;
            float dy = pa->array[i].y - pa->array[j].y;
            float dist2 = dx * dx + dy * dy;
            if (dist2 < 4.0f) { 
                pa->array[i].vx = -pa->array[i].vx;
                pa->array[i].vy = -pa->array[i].vy;
                pa->array[j].vx = -pa->array[j].vx;
                pa->array[j].vy = -pa->array[j].vy;
            }
        }
    }
}