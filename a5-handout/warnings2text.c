#include <stdio.h>
#include <stdlib.h>
#include "util.h"

int main(int argc, char** argv) {
  if (argc != 2) {
    fprintf(stderr, "Usage: %s FILE\n", argv[0]);
    exit(1);
  }

  const char* fname = argv[1];

  int32_t n;
  struct warning* ts = read_warnings(fname, &n);

  for (int i = 0; i < n; i++) {
    printf("%d %d\n", ts[i].s, ts[i].i);
  }
  free(ts);
}
