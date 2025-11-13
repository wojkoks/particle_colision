project(particle_collision C)
set(CMAKE_C_STANDARD 11)

# === OpenMP Support ===
find_package(OpenMP REQUIRED)

# === MPI Support ===
find_package(MPI REQUIRED)

# Target 1: OpenMP version
add_executable(particles src/main.c src/physics.c src/draw.c)

# Target 2: MPI version
add_executable(particles_mpi src/main_mpi.c src/physics_mpi.c src/draw.c)

find_package(PkgConfig REQUIRED)
pkg_check_modules(SDL2 REQUIRED sdl2)
include_directories(${SDL2_INCLUDE_DIRS})

# Link OpenMP executable
target_link_libraries(particles ${SDL2_LIBRARIES} OpenMP::OpenMP_C -lm)

# Link MPI executable
target_link_libraries(particles_mpi ${SDL2_LIBRARIES} MPI::MPI_C -lm)