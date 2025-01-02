#pragma once

#include <stdio.h>
#include <stdint.h>

// Return the number of seconds since some unspecified time. Timing
// can be done by calling this function multiple times and subtracting
// the return values.
double seconds(void);

// A vector in three-dimensional space. We use this for representing
// not just points, but also accelerations, velocities, etc.
struct vec3 {
  double x, y, z;
};

struct particle {
  struct vec3 pos;    // Positions.
  struct vec3 vel;    // Velocity.
  double mass;        // Mass.
};

struct warning {
  int s;  // Step
  int i;  // Index
};

// Compute Euclidean distance between two points in three-dimensional
// space.
double dist(struct vec3 a, struct vec3 b);

// Compute Euclidean distance between a point and (0,0,0) in three-dimensional
// space.
double dist_centre(struct vec3 v);

// Compute gravitational force of particle at position 'pj' with mass
// 'mj' on a particle at position 'pi'.
struct vec3 force(struct vec3 pi, struct vec3 pj, double mj);

// Read a particle file with name 'fname'. Returns an array of
// particles, and stores the number of particles in *n.
//
// This function does not report errors, but it *must* possible detect
// errors via assert().
struct particle* read_particles(const char *fname, int32_t *n);

// Write an array 'ps' of particles with the given size 'n' to a
// particle file with the name 'fname'.
//
// This function does not report errors, but it *must* possible detect
// errors via assert().
void write_particles(const char *fname, int32_t n, struct particle *ps);

// Read a warning file with name 'fname'. Returns an array of
// warnings, and stores the number of warnings in *n.
//
// This function does not report errors, but it *must* possible detect
// errors via assert().
struct warning* read_warnings(const char *fname, int32_t *n);

// Write an array 'ts' of warnings with the given size 'n' to a
// warning file with the name 'fname'.
//
// This function does not report errors, but it *must* possible detect
// errors via assert().
void write_warnings(const char *fname, int32_t n, struct warning *ts);

