#include <cstdio>
#include <string>
#include <stdlib.h>

#include "string_system.h"
#include "solver.h"

using std::string;

typedef struct _thread_args
{
  size_t thread_index;
  size_t thread_count;
  int argc;
  char** argv;
  string_system* my_system;
  pthread_barrier_t *sync;
} thread_args;

static void* solution_threaded (void* pa)
{
  thread_args *pargs = (thread_args*) pa;

  solve (pargs->thread_index, pargs->thread_count, pargs->argc, pargs->argv,
         pargs->my_system, pargs->sync);

  pthread_barrier_wait (pargs->sync);
  return 0;
}

static void process_system (
    size_t thread_count,
    int argc,
    char** argv,
    string_system* my_system
    )
{
  pthread_t index;
  thread_args *args = new thread_args [thread_count];
  pthread_barrier_t sync;
  pthread_barrier_init (&sync, 0, thread_count);

  for (size_t i = 0; i < thread_count; ++i)
    {
      args[i].argc = argc;
      args[i].argv = argv;
      args[i].my_system = my_system;
      args[i].sync = &sync;
      args[i].thread_count = thread_count;
      args[i].thread_index = i;
    }

  for (size_t i = 1; i < thread_count; ++i)
    {
      while (pthread_create (&index, NULL, solution_threaded, args + i))
        {
          fprintf (stderr, "Cannot create thread: %d\n", (int) i);
        }
    }
  solution_threaded (args + 0);

  pthread_barrier_destroy (&sync);
  delete[] args;
}

int main(int argc, char* argv[])
{
  int cpu_count = atol (argv[1]);

  string_system my_system;
  process_system (cpu_count, argc, argv, &my_system);
  my_system.print ();

  return 0;
}
