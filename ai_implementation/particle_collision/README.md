# Floating Particles Simulation

This project implements a particle simulation using C and SDL2, where a set of colorful particles move around the screen, bouncing off walls and each other. The simulation aims to run smoothly at a minimum of 60 frames per second.

## Project Structure

```
particle_collision
├── src
│   ├── main.c         # Entry point of the application
│   ├── physics.c      # Physics logic for particle simulation
│   ├── physics.h      # Header file for physics functions and structures
│   ├── draw.c         # Rendering logic for particles
│   └── draw.h         # Header file for rendering functions
├── CMakeLists.txt     # CMake configuration file
├── Makefile            # Makefile for building the project
└── README.md           # Project documentation
```

## Requirements

- C11 standard
- SDL2 library

## Building the Project

### Using CMake

1. Create a build directory:
   ```
   mkdir build && cd build
   ```

2. Run CMake:
   ```
   cmake ..
   ```

3. Build the project:
   ```
   cmake --build .
   ```

### Using Makefile

Simply run:
```
make
```

## Running the Simulation

After building the project, run the executable generated in the build directory. The simulation window will open, displaying the particles. Use the following controls:

- `Esc`: Close the application
- `Space`: Pause/Resume the simulation
- `R`: Restart with new random conditions
- `Up/Down`: Increase/Decrease the number of particles

## Functionality Overview

- **Particles**: Each particle has a random position, velocity, radius, and color. The simulation starts with at least 120 particles.
- **Collisions**: Particles bounce off walls and each other with elastic collisions.
- **Rendering**: Particles are drawn as filled circles on a black background.
- **Performance**: The simulation is designed to maintain a minimum of 60 FPS.

## License

This project is licensed under the MIT License. See the LICENSE file for more details.