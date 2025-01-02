#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include "util.h"

void print_particle(FILE* f, struct particle p) {
  fprintf(f, "%f %f %f %f %f %f %f\n",
          p.mass,
          p.pos.x, p.pos.y, p.pos.z,
          p.vel.x, p.vel.y, p.vel.z);
}

 // Feel free to fiddle with tolerance if necessary.
bool cmp(double x, double y) {
  double tol = 0.0000001;
  double rel_diff = fabs(fabs(x-y) / x);
  return rel_diff < tol;
}

int main(int argc, char** argv) {
  if (argc != 3) {
    fprintf(stderr, "Usage: %s FILE FILE\n", argv[0]);
    exit(1);
  }

  const char* fname_a = argv[1];
  const char* fname_b = argv[2];

  int32_t n, m;
  struct particle* ps_a = read_particles(fname_a, &n);
  struct particle* ps_b = read_particles(fname_b, &m);

  if (n != m) {
    fprintf(stderr, "Difference in particle count: %d != %d\n", n, m);
    return 1;
  }

  for (int i = 0; i < n; i++) {
    if (ps_a[i].mass != ps_b[i].mass ||
        !cmp(ps_a[i].pos.x, ps_b[i].pos.x) ||
        !cmp(ps_a[i].pos.y, ps_b[i].pos.y) ||
        !cmp(ps_a[i].pos.z, ps_b[i].pos.z) ||
        !cmp(ps_a[i].vel.x, ps_b[i].vel.x) ||
        !cmp(ps_a[i].vel.y, ps_b[i].vel.y) ||
        !cmp(ps_a[i].vel.z, ps_b[i].vel.z)) {
      fprintf(stderr, "Mismatch at particle %d\n", i);
      print_particle(stderr, ps_a[i]);
      print_particle(stderr, ps_b[i]);
      exit(1);
    }
  }
  free(ps_a);
  free(ps_b);
}
