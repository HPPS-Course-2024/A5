#include "util.h"
#include <assert.h>
#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>

// static const double WARNING_DISTANCE = 0.01;
static const double WARNING_DISTANCE = 0.01;

// Naive n-body simulation.
//
// *tc must be set to the number of warnings.
// *ts must point to an array of warnings with at least *tc elements.
void nbody(int n, struct particle* ps, int steps, int* tc,
           struct warning** ts) {
  *tc = 0; // Initialize total warning count
  *ts = NULL;

#pragma omp parallel
  {
    struct warning* thread_warnings = malloc(sizeof(struct warning));
    if (!thread_warnings) {
      fprintf(stderr, "malloc failed\n");
      exit(1);
    }
    int thread_warning_count    = 0;
    int thread_warning_capacity = 1;

#pragma omp for schedule(static)
    for (int s = 0; s < steps; s++) {
      for (int i = 0; i < n; i++) {
        double fx = 0, fy = 0, fz = 0;
        for (int j = 0; j < n; j++) {
          if (j == i) {
            continue;
          }
          struct vec3 a = force(ps[i].pos, ps[j].pos, ps[j].mass);

          fx += a.x;
          fy += a.y;
          fz += a.z;
        }
        ps[i].vel.x += fx;
        ps[i].vel.y += fy;
        ps[i].vel.z += fz;
      }
      for (int i = 0; i < n; i++) {
        ps[i].pos.x += ps[i].vel.x;
        ps[i].pos.y += ps[i].vel.y;
        ps[i].pos.z += ps[i].vel.z;
      }

      for (int i = 0; i < n; i++) {
        if (dist_centre(ps[i].pos) < WARNING_DISTANCE) {
          if (thread_warning_count >= thread_warning_capacity) {
            thread_warning_capacity *= 2;
            thread_warnings =
                realloc(thread_warnings,
                        thread_warning_capacity * sizeof(struct warning));
            if (!thread_warnings) {
              fprintf(stderr, "realloc failed\n");
              exit(1);
            }
          }
          printf("Warning: particle %d is too close to the centre\n", i);
          // thread_warnings[thread_warning_count - 1].s = s;
          // thread_warnings[thread_warning_count].i     = i;
          thread_warnings[thread_warning_count] = (struct warning){s, i};
          thread_warning_count++;
        }
      }
    }
#pragma omp critical
    {
      *ts = realloc(*ts, (*tc + thread_warning_count) * sizeof(struct warning));
      if (!*ts) {
        fprintf(stderr, "realloc failed\n");
        exit(1);
      }
      memcpy(*ts + *tc, thread_warnings,
             thread_warning_count * sizeof(struct warning));
      *tc += thread_warning_count;
    }
    free(thread_warnings);
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
