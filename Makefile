CFLAGS  = -std=c11 -O2 -Wall -Wextra -Wpedantic $(shell pkg-config --cflags sdl2)
LDFLAGS = $(shell pkg-config --libs sdl2) -lm

all: particles

particles: src/main.c src/physics.c src/draw.c
	$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS)

clean:
	rm -f particles
