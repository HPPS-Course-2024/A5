#include "util.h"
#include <assert.h>
#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
// static const double WARNING_DISTANCE = 0.01;
static const double WARNING_DISTANCE = 0.01;

// Naive n-body simulation.
//
// *tc must be set to the number of warnings.
// *ts must point to an array of warnings with at least *tc elements.
void nbody(int n, struct particle* ps, int steps, int* tc,
           struct warning** ts) {
  // Start with an empty global warning array
  *tc = 0;
  *ts = NULL;
  printf("n: %d\n", n);
  // Cannot parallelize this loop as, you don’t parallelize the outer steps loop
  // in an N-body simulation because each step depends on the updated positions
  // and velocities from the previous step. That means step s+1 can’t start
  // until step s has finished. So the outer loop is usually left sequential,
  // and you parallelize the inner particle loops for forces, positions, etc.
  for (int s = 0; s < steps; s++) {
// 1) Parallel loop for forces
#pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < n; i++) {
      double fx = 0, fy = 0, fz = 0;
      for (int j = 0; j < n; j++) {
        if (j != i) {
          struct vec3 a = force(ps[i].pos, ps[j].pos, ps[j].mass);
          fx += a.x;
          fy += a.y;
          fz += a.z;
        }
      }
      // Update velocity
      ps[i].vel.x += fx;
      ps[i].vel.y += fy;
      ps[i].vel.z += fz;
    }

// 2) Parallel loop for positions
#pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < n; i++) {
      ps[i].pos.x += ps[i].vel.x;
      ps[i].pos.y += ps[i].vel.y;
      ps[i].pos.z += ps[i].vel.z;
    }

// 3) Collect warnings privately per thread, then merge
#pragma omp parallel
    {
      // Use a local dynamic array (or something similar) in each thread
      struct warning* local_warnings = NULL;
      int             local_cap      = 0;
      int             local_count    = 0;

// Each thread loops over all particles (can use #pragma omp for again)
#pragma omp for schedule(dynamic)
      for (int i = 0; i < n; i++) {
        double distance = dist_centre(ps[i].pos);
        if (distance < WARNING_DISTANCE) {
          printf("<WARNING> Particle %d is too close to the centre!\n", i);
          if (local_count >= local_cap) {
            local_cap = (local_cap == 0) ? 1 : local_cap * 2;
            local_warnings =
                realloc(local_warnings, local_cap * sizeof(struct warning));
            if (!local_warnings) {
              fprintf(stderr, "realloc failed\n");
              exit(1);
            }
          }
          // Add warning
          local_warnings[local_count++] = (struct warning){s, i};
        }
      }

// Merge local warnings into the global array
#pragma omp critical
      {
        if (local_count > 0) {
          int old_tc = *tc;
          *tc += local_count;
          *ts = realloc(*ts, (*tc) * sizeof(struct warning));
          if (!*ts) {
            fprintf(stderr, "realloc failed\n");
            exit(1);
          }
          memcpy((*ts) + old_tc, local_warnings,
                 local_count * sizeof(struct warning));
        }
      }

      free(local_warnings);
    } // end parallel warning collection
  }
}

int main(int argc, char** argv) {
  int steps = 1;
  if (argc < 4) {
    printf("Usage: \n");
    printf("%s <input> <particle output> <warnings output> [steps]\n", argv[0]);
    return 1;
  } else if (argc > 4) {
    steps = atoi(argv[4]);
  }

  int32_t          n;
  struct particle* ps = read_particles(argv[1], &n);

  int             tc = 0;    // int to store size of warning array
  struct warning* ts = NULL; // array to store warnings

  double bef = seconds();
  nbody(n, ps, steps, &tc, &ts);
  double aft = seconds();
  printf("%f\n", aft - bef);
  write_particles(argv[2], n, ps);
  write_warnings(argv[3], tc, ts);

  free(ts);
  free(ps);
}
