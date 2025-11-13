#include <SDL2/SDL.h>
#include <math.h>

void draw_filled_circle(SDL_Renderer *r, int cx, int cy, int radius)
{
    for (int dy = -radius; dy <= radius; ++dy)
    {
        int dx = (int)floorf(sqrtf((float)radius * radius - (float)dy * dy));
        SDL_RenderDrawLine(r, cx - dx, cy + dy, cx + dx, cy + dy);
    }
}

void draw_grid(SDL_Renderer *r, int width, int height, int cols, int rows)
{
    // Draw vertical lines
    int col_width = width / cols;
    SDL_SetRenderDrawColor(r, 100, 100, 100, 255);  // Gray
    for (int col = 1; col < cols; col++)
    {
        SDL_RenderDrawLine(r, col * col_width, 0, col * col_width, height);
    }
    
    // Draw horizontal lines
    int row_height = height / rows;
    for (int row = 1; row < rows; row++)
    {
        SDL_RenderDrawLine(r, 0, row * row_height, width, row * row_height);
    }
}