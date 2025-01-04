#include "util.h"
#include <assert.h>
#include <float.h>
#include <math.h>
#include <omp.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define NUM_OCTANTS 8

int octant(struct vec3 corner, double l, const struct particle* p) {
  if (p->pos.x >= corner.x + l / 2) {
    if (p->pos.y >= corner.y + l / 2) {
      if (p->pos.z >= corner.z + l / 2) {
        return 0;
      } else {
        return 1;
      }
    } else {
      if (p->pos.z >= corner.z + l / 2) {
        return 2;
      } else {
        return 3;
      }
    }
  } else {
    if (p->pos.y >= corner.y + l / 2) {
      if (p->pos.z >= corner.z + l / 2) {
        return 4;
      } else {
        return 5;
      }
    } else {
      if (p->pos.z >= corner.z + l / 2) {
        return 6;
      } else {
        return 7;
      }
    }
  }
}

void octant_offset(int j, double* ox, double* oy, double* oz) {
  switch (j) {
    case 0:
      *ox = 0.5;
      *oy = 0.5;
      *oz = 0.5;
      break;
    case 1:
      *ox = 0.5;
      *oy = 0.5;
      *oz = 0.0;
      break;
    case 2:
      *ox = 0.5;
      *oy = 0.0;
      *oz = 0.5;
      break;
    case 3:
      *ox = 0.5;
      *oy = 0.0;
      *oz = 0.0;
      break;
    case 4:
      *ox = 0.0;
      *oy = 0.5;
      *oz = 0.5;
      break;
    case 5:
      *ox = 0.0;
      *oy = 0.5;
      *oz = 0.0;
      break;
    case 6:
      *ox = 0.0;
      *oy = 0.0;
      *oz = 0.5;
      break;
    case 7:
      *ox = 0.0;
      *oy = 0.0;
      *oz = 0.0;
      break;
  }
}

struct bh_node {
  bool internal; // False when external.

  struct vec3 corner;
  double      l; // Edge length of space.

  // Fields for external nodes.
  int particle; // Index of particle in particle array; -1 if none.

  // Fields for internal nodes. Only have sensible values when
  // 'internal' is true.
  struct vec3     com;  // Center of mass.
  double          mass; // Total mass.
  struct bh_node* children[8];
};

void set_vec3_fields(struct vec3* v, double x, double y, double z) {
  v->x = x;
  v->y = y;
  v->z = z;
}

void update_center_of_mass(struct bh_node* bh) {
  assert(bh->internal);
  double      total_mass = 0.0;
  struct vec3 com        = {0, 0, 0};

  for (int i = 0; i < NUM_OCTANTS; i++) {
    struct bh_node* child = bh->children[i];
    if (!child) {
      fprintf(stderr, "Child is NULL\n");
      continue;
    }
    double child_mass = child->mass;
    total_mass += child_mass;
    com.x += child->com.x * child_mass;
    com.y += child->com.y * child_mass;
    com.z += child->com.z * child_mass;
  }
  bh->mass = total_mass;

  if (total_mass > 0.0) {
    com.x /= total_mass;
    com.y /= total_mass;
    com.z /= total_mass;
    bh->com = com;
  }
}

void bh_mk_internal(struct bh_node* bh) {
  // Must not already be internal.
  assert(!bh->internal);

  bh->particle = -1;
  bh->internal = true;
  bh->mass     = 0;
  set_vec3_fields(&bh->com, 0, 0, 0);

  for (int i = 0; i < NUM_OCTANTS; i++) {
    bh->children[i] = malloc(sizeof(struct bh_node));
    if (!bh->children[i]) {
      fprintf(stderr, "malloc failed for children\n");
      exit(EXIT_FAILURE);
    }
    struct bh_node* child = bh->children[i];
    child->internal       = false;
    child->particle       = -1;
    child->mass           = 0;
    set_vec3_fields(&child->com, 0, 0, 0);

    double ox, oy, oz;
    octant_offset(i, &ox, &oy, &oz);

    double corner_x = bh->corner.x + bh->l * ox;
    double corner_y = bh->corner.y + bh->l * oy;
    double corner_z = bh->corner.z + bh->l * oz;
    set_vec3_fields(&child->corner, corner_x, corner_y, corner_z);

    child->l = bh->l / 2.0;
  }
}

void bh_insert(struct bh_node* bh, struct particle* ps, int p) {
  if (bh->internal) {
    int oct = octant(bh->corner, bh->l, ps + p);
    bh_insert(bh->children[oct], ps, p);
    update_center_of_mass(bh);
  } else {
    if (bh->particle == -1) {
      bh->particle = p;
      bh->mass     = ps[p].mass;
      bh->com      = ps[p].pos;
    } else {
      int old_p_idx = bh->particle;
      bh_mk_internal(bh);

      int oct_old = octant(bh->corner, bh->l, ps + old_p_idx);
      bh_insert(bh->children[oct_old], ps, old_p_idx);

      int oct_new = octant(bh->corner, bh->l, ps + p);
      bh_insert(bh->children[oct_new], ps, p);

      update_center_of_mass(bh);
    }
  }
}

void bh_free(struct bh_node* bh) {
  if (bh->internal) {
    for (int i = 0; i < NUM_OCTANTS; i++) {
      bh_free(bh->children[i]);
    }
  }
  free(bh);
}

void bh_accel(double theta, struct bh_node* node, struct particle* ps, int p,
              struct vec3* a) {
  if (!node)
    return; // Safety check, though node should never be NULL if used correctly.

  if (node->internal) {
    double d = dist(node->com, (ps + p)->pos);
    if (d > 0.0 && node->l / d < theta) {
      struct vec3 f = force((ps + p)->pos, node->com, node->mass);
      a->x += f.x;
      a->y += f.y;
      a->z += f.z;
    } else {
      for (int i = 0; i < NUM_OCTANTS; i++) {
        bh_accel(theta, node->children[i], ps, p, a);
      }
    }
  } else if (node->particle != -1 && node->particle != p) {
    struct particle* p_node = ps + node->particle;
    struct vec3      f      = force((ps + p)->pos, p_node->pos, p_node->mass);
    a->x += f.x;
    a->y += f.y;
    a->z += f.z;
  }
}

struct bh_node* bh_new(double min_coord, double max_coord) {
  struct bh_node* bh = malloc(sizeof(struct bh_node));
  if (!bh) {
    fprintf(stderr, "malloc failed in bh_new()\n");
    exit(EXIT_FAILURE);
  }
  bh->corner.x = min_coord;
  bh->corner.y = min_coord;
  bh->corner.z = min_coord;
  bh->l        = max_coord - min_coord;
  bh->internal = false;
  bh->particle = -1;
  bh->mass     = 0.0;
  set_vec3_fields(&bh->com, 0, 0, 0);

  // Initialize children pointers to NULL to avoid any uninitialized usage
  for (int i = 0; i < NUM_OCTANTS; i++) {
    bh->children[i] = NULL;
  }
  return bh;
}

static const double WARNING_DISTANCE = 0.01;

void nbody(int n, struct particle* ps, int steps, int* tc, struct warning** ts,
           double theta) {
  for (int s = 0; s < steps; s++) {
    double max_coord = -DBL_MAX;
    double min_coord = DBL_MAX;

#pragma omp parallel for reduction(max : max_coord) reduction(min : min_coord)
    for (int i = 0; i < n; i++) {
      double px = ps[i].pos.x;
      double py = ps[i].pos.y;
      double pz = ps[i].pos.z;

      if (px < min_coord)
        min_coord = px;
      if (px > max_coord)
        max_coord = px;

      if (py < min_coord)
        min_coord = py;
      if (py > max_coord)
        max_coord = py;

      if (pz < min_coord)
        min_coord = pz;
      if (pz > max_coord)
        max_coord = pz;
    }

    struct bh_node* root = bh_new(min_coord, max_coord);

#pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < n; i++) {
#pragma omp critical
      { bh_insert(root, ps, i); }
    }

#pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < n; i++) {
      struct vec3 a = {0, 0, 0};
      bh_accel(theta, root, ps, i, &a);
      ps[i].vel.x += a.x;
      ps[i].vel.y += a.y;
      ps[i].vel.z += a.z;
    }

#pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < n; i++) {
      ps[i].pos.x += ps[i].vel.x;
      ps[i].pos.y += ps[i].vel.y;
      ps[i].pos.z += ps[i].vel.z;
    }

#pragma omp parallel
    {
      int             local_cap   = 1;
      int             local_count = 0;
      struct warning* local_warnings =
          malloc(local_cap * sizeof(struct warning));
      if (!local_warnings) {
        fprintf(stderr, "malloc failed\n");
        exit(EXIT_FAILURE);
      }

#pragma omp for schedule(dynamic)
      for (int i = 0; i < n; i++) {
        double distance = dist_centre(ps[i].pos);
        if (distance < WARNING_DISTANCE) {
          if (local_count >= local_cap) {
            local_cap *= 2;
            local_warnings =
                realloc(local_warnings, local_cap * sizeof(struct warning));
            if (!local_warnings) {
              fprintf(stderr, "realloc failed\n");
              exit(EXIT_FAILURE);
            }
          }
          printf("<WARNING> Particle %d is too close to the centre!\n", i);
          local_warnings[local_count++] = (struct warning){s, i};
        }
      }

#pragma omp critical
      {
        if (local_count > 0) {
          int old_tc = *tc;
          *tc += local_count;
          *ts = realloc(*ts, (*tc) * sizeof(struct warning));
          if (!*ts) {
            fprintf(stderr, "realloc failed\n");
            exit(EXIT_FAILURE);
          }
          memcpy((*ts) + old_tc, local_warnings,
                 local_count * sizeof(struct warning));
        }
      }
      free(local_warnings);
    }

    bh_free(root);
  }
}

int main(int argc, char** argv) {
  int    steps = 1;
  double theta = 0.5;
  if (argc < 4) {
    printf("Usage: \n");
    printf("%s <input> <particle output> <warnings output> [steps] [theta]\n",
           argv[0]);
    return 1;
  }
  if (argc > 4) {
    steps = atoi(argv[4]);
  }
  if (argc > 5) {
    theta = atof(argv[5]);
  }

  int32_t          n;
  struct particle* ps = read_particles(argv[1], &n);

  int             tc = 0;    // int to store size of warning array
  struct warning* ts = NULL; // array to store warnings

  double bef = seconds();
  nbody(n, ps, steps, &tc, &ts, theta);
  double aft = seconds();
  printf("%f\n", aft - bef);

  write_particles(argv[2], n, ps);
  write_warnings(argv[3], tc, ts);

  free(ts);
  free(ps);

  return 0;
}
