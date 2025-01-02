#include "util.h"
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

double seconds(void) {
  struct timeval tv;
  gettimeofday(&tv, NULL); // The NULL is for timezone information.
  return tv.tv_sec + tv.tv_usec / 1000000.0;
}

double dist(struct vec3 a, struct vec3 b) {
  double dx = a.x - b.x;
  double dy = a.y - b.y;
  double dz = a.z - b.z;

  return sqrt(dx * dx + dy * dy + dz * dz);
}
//(x1-x2)^2+(y1-y2)^2+(z1-z2)^2
double dist_centre(struct vec3 v) {
  return sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
}

struct vec3 force(struct vec3 pi, struct vec3 pj, double mj) {
  static const double epsilon = 1;

  double d      = dist(pi, pj);
  double d_sqr  = d * d;
  double inv_d  = 1.0 / sqrt(d_sqr + epsilon * epsilon);
  double inv_d3 = inv_d * inv_d * inv_d;

  struct vec3 f;
  f.x = (pj.x - pi.x) * mj * inv_d3;
  f.y = (pj.y - pi.y) * mj * inv_d3;
  f.z = (pj.z - pi.z) * mj * inv_d3;
  return f;
}

struct particle* read_particles(const char* fname, int32_t* n) {
  FILE* f = fopen(fname, "r");
  assert(f != NULL);

  assert(fread(n, sizeof(int32_t), 1, f) == 1);

  struct particle* ps = malloc(*n * sizeof(struct particle));
  assert(ps != NULL);

  for (int i = 0; i < *n; i++) {
    double x[7];
    assert(fread(x, sizeof(double), 7, f) == 7);
    ps[i].mass  = x[0];
    ps[i].pos.x = x[1];
    ps[i].pos.y = x[2];
    ps[i].pos.z = x[3];
    ps[i].vel.x = x[4];
    ps[i].vel.y = x[5];
    ps[i].vel.z = x[6];
  }

  fclose(f);

  return ps;
}

void write_particles(const char* fname, int32_t n, struct particle* ps) {
  FILE* f = fopen(fname, "w");
  assert(f != NULL);

  assert(fwrite(&n, sizeof(int32_t), 1, f) == 1);

  for (int i = 0; i < n; i++) {
    assert(fwrite(&ps[i].mass, sizeof(double), 1, f) == 1);
    assert(fwrite(&ps[i].pos.x, sizeof(double), 1, f) == 1);
    assert(fwrite(&ps[i].pos.y, sizeof(double), 1, f) == 1);
    assert(fwrite(&ps[i].pos.z, sizeof(double), 1, f) == 1);
    assert(fwrite(&ps[i].vel.x, sizeof(double), 1, f) == 1);
    assert(fwrite(&ps[i].vel.y, sizeof(double), 1, f) == 1);
    assert(fwrite(&ps[i].vel.z, sizeof(double), 1, f) == 1);
  }

  fclose(f);
}

struct warning* read_warnings(const char* fname, int32_t* n) {
  FILE* f = fopen(fname, "r");
  assert(f != NULL);

  assert(fread(n, sizeof(int), 1, f) == 1);

  struct warning* ts = malloc(*n * sizeof(struct warning));
  assert(ts != NULL);

  for (int i = 0; i < *n; i++) {
    int x[2];
    assert(fread(x, sizeof(int), 2, f) == 2);
    ts[i].s = x[0];
    ts[i].i = x[1];
  }

  fclose(f);

  return ts;
}

void write_warnings(const char* fname, int32_t n, struct warning* ts) {
  FILE* f = fopen(fname, "w");
  assert(f != NULL);

  assert(fwrite(&n, sizeof(int), 1, f) == 1);

  for (int i = 0; i < n; i++) {
    assert(fwrite(&ts[i].s, sizeof(int), 1, f) == 1);
    assert(fwrite(&ts[i].i, sizeof(int), 1, f) == 1);
  }

  fclose(f);
}
