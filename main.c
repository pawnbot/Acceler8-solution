#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sched.h>
#include <pthread.h>


#include "timer.h"
#include "array_op.h"
#include "kadane.h"
#include "processor_op.h"

// 128 x 128 matrix
#define SIZE1 65536
// L3 cache size
#define SIZE2 25165824
// Debug output
// #define LOCAL_DEBUG

typedef struct _ARGS_MAIN /* main */
{
  int n;                          // rows in matrix
  int m;                          // columns in matrix
  int total_threads;              // number of threads
  int* ps;                        // (n x m) sums matrix
  int* workspace;                 // (max(MIN_LEN, M) + 2 * M) workspace
  answer *result;                 // result
  
  // set cpu
  int start_thread;
  int active_threads;
  // pseudo random  
  int A; 
  int B; 
  int M; 
  int SEED0;
} ARGS_MAIN;

static __inline__ int
max(int a, int b) {
  return (a > b ? a: b);
}

static __inline__ int
min(int a, int b) {
  return (a < b ? a: b);
}

static void * main_threaded(void *pargs)
{
  ARGS_MAIN *pa = (ARGS_MAIN*) pargs;

  kadane_alg(pa->n, pa->m, pa->total_threads,
             pa->start_thread, pa->active_threads,
             pa->ps, pa->workspace, pa->result, 
             pa->A, pa->B, pa->M, pa->SEED0);
  
  return 0;
}

int main(int argc, char *argv[]) {
  //////////////////////////////////////////////////////////////////////////////
  int active_threads;           // total number of active threads
  int case_threads;             // threads on case
  //////////////////////////////////////////////////////////////////////////////  
  int total_tests;              // Total number of tests
  int row_size, col_size;       // number of rows and columns in matrix
  
  int iter;                     // loop
  int a, b, m, seed0;           // pseudo random
  int ret_val;                  // return value
  //////////////////////////////////////////////////////////////////////////////
  // Input parallelization
  ARGS_MAIN *args;
  pthread_t *threads;
  answer *results;
  int *ps;
  size_t begin, k, total_mem;
  //////////////////////////////////////////////////////////////////////////////
  FILE *input, *output;         // input, output
  long int timer;

#ifdef LOCAL_DEBUG  
  //////////////////////////////////////////////////////////////////////////////
  // Timer started
  timer_start ();
  timer = get_full_time();
#endif /* LOCAL_DEBUG */    
    
  get_proc_info(&active_threads, 0, 0);
  printf("Logical CPUs detected: %d\n", 
         active_threads);
  
  if ((input = fopen (argv[1], "r")) == NULL) {
    printf("Cannot open input file!\n");  
    return 1;
  }

  ret_val = fscanf(input, "%d", &total_tests);
  while (!(results = (answer*) malloc(total_tests * sizeof(answer))))
    printf("Waiting for memory for results...\n");
  while (!(args = (ARGS_MAIN*) malloc(active_threads * sizeof(ARGS_MAIN)))) 
    printf("Waiting for memory for args...\n");
  while (!(threads = (pthread_t*) malloc(active_threads * sizeof(pthread_t)))) 
    printf("Waiting for memory for threads...\n");
  
  memset(threads, 0, active_threads * sizeof(pthread_t));
  begin = 0;
  for (iter = 0; iter < total_tests; ++iter) {
    ret_val = fscanf(input, "%d %d", &row_size, &col_size);
    ret_val = fscanf(input, "%d %d %d %d", &seed0, &a, &b, &m);
    // a_{i j} = const
    if ((a * seed0 + b) % m == seed0) {
      results[iter].left_top_row = 0;
      results[iter].left_top_col = 0;
      results[iter].right_bot_row = 0;
      results[iter].right_bot_col = 0;
      results[iter].sum = 0;
      results[iter].square = 1;
      continue;
    }

    total_mem = ((size_t) row_size * col_size + max(MIN_LEN, m) + 2 * m) * sizeof(int);
    case_threads = min(active_threads, (row_size + BLOCK_ROW - 1) / BLOCK_ROW);
    
    if (total_mem < SIZE1) {
      case_threads = min(case_threads, 1); 
    } else {
      if (total_mem < SIZE2 && total_tests > 7)
        case_threads = min(case_threads, 10);
    }
    
    for (k = 0; k < case_threads; ++k) {
      if (threads[(begin + k) % active_threads] != 0)
        pthread_join(threads[(begin + k) % active_threads], NULL);
      threads[(begin + k) % active_threads] = 0;
    }
    k = case_threads;
    while (!(ps = (int*) malloc(total_mem))) {
      if (threads[(begin + k) % active_threads] != 0)
        pthread_join(threads[(begin + k) % active_threads], NULL);
      threads[(begin + k) % active_threads] = 0;
      ++k;
    }

    args[begin].n = row_size;
    args[begin].m = col_size;
    args[begin].total_threads = case_threads;
    args[begin].SEED0 = seed0;
    args[begin].A = a;
    args[begin].B = b;
    args[begin].M = m;
    args[begin].result = results + iter;
    args[begin].ps = ps;
    args[begin].workspace = ps + row_size * col_size;
    args[begin].start_thread = begin;
    args[begin].active_threads = active_threads;
    
    while (pthread_create (threads + begin, NULL, 
                           main_threaded,
                           args + begin))
    {
      if (threads[(begin + k) % active_threads] != 0)
        pthread_join(threads[(begin + k) % active_threads], NULL);
      threads[(begin + k) % active_threads] = 0;
      ++k;
    }
    begin = (begin + case_threads) % active_threads;

#ifdef LOCAL_DEBUG
    printf("Case %d, Threads: %d, Started in: %ld.%ld (sec)\n", 
           iter + 1, case_threads, (get_full_time() - timer)/100, 
           (get_full_time() - timer) % 100);
#endif /* LOCAL_DEBUG */
  }
  
  for (k = 0; k < active_threads; ++k) {
    if (threads[(begin + k) % active_threads] != 0)
      pthread_join(threads[(begin + k) % active_threads], NULL);
    threads[(begin + k) % active_threads] = 0;
  }
  free(args);
  free(threads);
  fclose(input);

  if ((output = fopen (argv[2], "w")) == NULL) {
    printf("Cannot open output file!\n");  
    return 2;
  }
  
  for (iter = 0; iter < total_tests; ++iter) {
    ret_val = fprintf(output, "Case #%d: %d %d %d %d %d %d\n", iter + 1,
                      results[iter].left_top_row, results[iter].left_top_col,
                      results[iter].right_bot_row, results[iter].right_bot_col,
                      results[iter].sum, results[iter].square);
  }
  
  free(results);
  fclose(output);
  
#ifdef LOCAL_DEBUG    
  //////////////////////////////////////////////////////////////////////////////
  // Print full time
  timer = get_full_time() - timer;
  printf("Total time: %ld.%ld (sec)\n", timer / 100, timer % 100);
#endif /* LOCAL_DEBUG */
    
  return 0;
}

