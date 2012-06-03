#ifndef KADANE_H
#define KADANE_H

typedef struct _ARGS /* kadane_alg */
{
  int n;                   // rows in matrix
  int m;                   // columns in matrix
  int total_threads;       // number of threads
  int thread_number;       // thread number 
  int* ps;                 // (n x m) sums matrix
  int* workspace;          // (2 x BLOCK_ROW x BLOCK_COL + BLOCK_ROW +
                           // BLOCK_ROW x BLOCK_ROW) workspace
  
  answer result;           // result
  pthread_barrier_t* sync; // synchronize barrier
  
  int MAX_ROW;
  int CPU_ID;
  // cycle constants
  int *cycle, *precycle;
  int len, prelen;    
} ARGS;

/*******************************************************************************
* Kadane algorithm
*******************************************************************************/
int kadane_alg (
  int n,                          // rows in matrix
  int m,                          // columns in matrix
  int total_threads,              // number of threads
  int start_thread,               // Global begin thread
  int active_threads,             // Global total threads
  int* ps,                        // (n x m) sums matrix
  int* workspace,                 // (max(MIN_LEN, M) + 2 * M) workspace
  answer *result,                 // result
  int A, int B, int M, int SEED0  // pseudo random
  );

#endif /* KADANE_H */
