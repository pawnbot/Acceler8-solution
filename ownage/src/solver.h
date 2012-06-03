#ifndef SOLVER_H
#define SOLVER_H
#include "string_system.h"

void solve (
  size_t thread_index,
  size_t thread_count,
  int argc,
  char* argv[],
  string_system* my_system,
  pthread_barrier_t *sync
  );

#endif // SOLVER_H
