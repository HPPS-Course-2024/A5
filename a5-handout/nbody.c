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
  int warning_capacity = 1; // Start with capacity 1
  *ts                  = malloc(warning_capacity * sizeof(struct warning));
  if (!*ts) {
    fprintf(stderr, "malloc failed\n");
    exit(1);
  }

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
      double distance = dist_centre(ps[i].pos);
      if (distance < WARNING_DISTANCE) {
        if (*tc >= warning_capacity) {
          warning_capacity *= 2;
          *ts = realloc(*ts, warning_capacity * sizeof(struct warning));
          if (!*ts) {
            fprintf(stderr, "realloc failed\n");
            exit(1);
          }
        }
        (*ts)[*tc].s = s;
        (*ts)[*tc].i = i;
        (*tc)++;
      }
    }
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
  printf("%fpenis\n", aft - bef);
  write_particles(argv[2], n, ps);
  write_warnings(argv[3], tc, ts);

  free(ts);
  free(ps);
}
