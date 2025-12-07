#include <SDL2/SDL.h>
#include <math.h>
#include <string.h>

typedef struct
{
    char ch;
    uint8_t rows[5]; // 3 bits used
} Glyph3x5;

static const Glyph3x5 FONT_3X5[] = {
    {'0', {7, 5, 5, 5, 7}},
    {'1', {2, 6, 2, 2, 7}},
    {'2', {7, 1, 7, 4, 7}},
    {'3', {7, 1, 7, 1, 7}},
    {'4', {5, 5, 7, 1, 1}},
    {'5', {7, 4, 7, 1, 7}},
    {'6', {7, 4, 7, 5, 7}},
    {'7', {7, 1, 1, 1, 1}},
    {'8', {7, 5, 7, 5, 7}},
    {'9', {7, 5, 7, 1, 7}},
    {'F', {7, 4, 7, 4, 4}},
    {'P', {7, 5, 7, 4, 4}},
    {'S', {7, 4, 7, 1, 7}},
    {'.', {0, 0, 0, 0, 1}},
    {' ', {0, 0, 0, 0, 0}},
};

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

static const uint8_t *find_glyph(char c)
{
    size_t count = sizeof(FONT_3X5) / sizeof(FONT_3X5[0]);
    for (size_t i = 0; i < count; i++)
    {
        if (FONT_3X5[i].ch == c)
            return FONT_3X5[i].rows;
    }
    return FONT_3X5[count - 1].rows; // space
}

static void draw_glyph(SDL_Renderer *r, int x, int y, char c, int scale, SDL_Color color)
{
    const uint8_t *rows = find_glyph(c);
    SDL_SetRenderDrawColor(r, color.r, color.g, color.b, color.a);
    for (int row = 0; row < 5; row++)
    {
        for (int col = 0; col < 3; col++)
        {
            if (rows[row] & (1 << (2 - col)))
            {
                SDL_Rect px = {x + col * scale, y + row * scale, scale, scale};
                SDL_RenderFillRect(r, &px);
            }
        }
    }
}

void draw_text(SDL_Renderer *r, int x, int y, const char *text, int scale, SDL_Color color)
{
    int cursor_x = x;
    size_t len = strlen(text);
    for (size_t i = 0; i < len; i++)
    {
        draw_glyph(r, cursor_x, y, text[i], scale, color);
        cursor_x += (3 * scale) + scale; // glyph width + spacing
    }
}
