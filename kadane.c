#include <stdlib.h>
#include <stdio.h>
#include <pthread.h>
#include <stdint.h>
#include <string.h>

#include "array_op.h"
#include "kadane.h"
#include "processor_op.h"

static __inline__ int
max(int a, int b) {
  return (a > b ? a: b);
}

static __inline__ int
min(int a, int b) {
  return (a < b ? a: b);
}

static __inline__ void
put_max(int* a, int b) {
  if (*a < b)
    *a = b;
}

static __inline__ void
put_min(int* a, int b) {
  if (*a > b)
    *a = b;
}

// swap integer A with B
static __inline__ void
swap (int* A, int* B) {
  int buff;
  buff = *A;
  *A = *B;
  *B = buff;
}

/*******************************************************************************
* Kadane algorithm
*******************************************************************************/
static void kadane_thread (
  int n,                   // rows in matrix
  int m,                   // columns in matrix
  int MAX_ROW,             // max top row of submatrix
  int total_threads,       // number of threads
  int thread_number,       // thread number  
  int *ps,                 // (n x m) sums matrix
  int* workspace,          // (2 x BLOCK_ROW x BLOCK_COL + BLOCK_ROW +
                           // BLOCK_ROW x BLOCK_ROW) workspace
  answer *result           // result
  )
{
  const int MIN_SUM = -20001;
  //////////////////////////////////////////////////////////////////////////////
  // loop variables
  size_t ltr, rbr, rbc;
  size_t rs, re, ks, ke, cs, ce, cl, iter, i, j;
  // resulting variables
  int sum;
  int res_rbc, res_ltc;
  int max_sum;
  int min_sum, min_ltc;
  answer res;
  // cache
  int *pps;         // pointer on ps
  int *top;         // top block of submatrix
  int *bot;         // bot block of submatrix
  int *pbot;        // pointer on bot
  int *ptop;        // pointer on top
  int *maxs;        // maximal sum array (BLOCK_ROW)
  int *mins;        // mininal sum array (BLOCK_ROW x BLOCK_ROW)
  int *pmins;       // pointer on msums
  
  // Initialization
  mins = workspace;
  maxs = mins + BLOCK_ROW * BLOCK_ROW;
  top = maxs + BLOCK_ROW;
  bot = top + BLOCK_ROW * BLOCK_COL;
  
  max_sum = MIN_SUM;
  res.sum = MIN_SUM;
  res.left_top_row = 0;
  res.left_top_col = 0;
  res.right_bot_row = 0;
  res.right_bot_col = 0;  
  res_rbc = 0;
  res_ltc = 0;
  for (i = 0; i < BLOCK_ROW; ++i)
    maxs[i] = MIN_SUM;
  
  // zero row
  ltr = 0;
  for (rbr = ltr + thread_number; rbr < n; rbr += total_threads) {
    // Init
    min_sum = 0;
    min_ltc = 0;
    pps = ps + rbr*m;
    sum = 0;
    
    for (rbc = 0; rbc < m; ++rbc) {
      sum = pps[rbc];
      if (sum - min_sum > max_sum) {
        max_sum = sum - min_sum;
        res_rbc = rbc;
        res_ltc = min_ltc;
      }
      if (sum < min_sum) {
        min_sum = sum;
        min_ltc = rbc + 1;
      }    
    }
    
    // Check_res
    if (max_sum > res.sum) {
      res.sum = max_sum;
      res.left_top_row = ltr; 
      res.left_top_col = res_ltc;
      res.right_bot_row = rbr; 
      res.right_bot_col = res_rbc;
    }
  }
  // zero row computed
  // main loop
  for (rs = 1; rs < MAX_ROW; rs += BLOCK_ROW) 
  {
    re = (rs + BLOCK_ROW < MAX_ROW ? rs + BLOCK_ROW : MAX_ROW);
    for (ks = rs + thread_number * BLOCK_ROW; 
         ks < n; ks += total_threads * BLOCK_ROW) 
    {
      ke = (ks + BLOCK_ROW < n ? ks + BLOCK_ROW : n);
      // set min sums to zero
      memset(mins, 0, BLOCK_ROW * BLOCK_ROW * sizeof(int));
      
      for (cs = 0; cs < m; cs += BLOCK_COL) {
        ce = (cs + BLOCK_COL < m ? cs + BLOCK_COL : m);
        cl = ce - cs;
        // copy and transpose bot block of submatrix
        pps  = ps + ks * m + cs;
        for (i = 0; i < ke - ks; ++i, pps += m) {
          pbot = bot + i;
          for (j = 0; j < cl; ++j)
            pbot[BLOCK_ROW * j] = pps[j];
        } 
        for (ltr = rs; ltr < re; ++ltr) {
          ptop = top + (ltr - rs) * BLOCK_COL;
          memcpy(ptop, ps + (ltr - 1) * m + cs, cl * sizeof(int));
          pmins = mins + (ltr - rs) * BLOCK_ROW;
          
          iter = ke - max(ltr, ks);
          pbot = bot + (max(ltr, ks) - ks);
          for (rbc = 0; rbc < cl; ++rbc, pbot += BLOCK_ROW) {
            sum = ptop[rbc];
            ////////////////vectorized main loop////////////////////////////////
#pragma simd
#pragma vector always
            for (rbr = 0; rbr < iter; ++rbr) {
              put_max(&maxs[rbr], pbot[rbr]  - pmins[rbr] - sum);
              put_min(&pmins[rbr], pbot[rbr] - sum); 
            }
            //////////////////////////////////////////////////////////////////// 
          }
          for (i = 0; i < BLOCK_ROW; ++i) {
            put_max(&max_sum, maxs[i]);
            if (max_sum > res.sum) {
              res.sum = max_sum;
              res.left_top_row = ltr; 
              res.right_bot_row = max(ltr, ks) + i; 
            }
          }
        }
      }
    }
  }
  
  // Calculating column numbers
  if (res.left_top_row > 0) {
    ltr = res.left_top_row;
    rbr = res.right_bot_row;
    // Init
    sum = 0;
    min_sum = 0;
    min_ltc = 0;
    ptop = ps + (ltr - 1)*m;
    pps = ps + rbr*m;

    for (rbc = 0; rbc < m; ++rbc) {
      sum = pps[rbc] - ptop[rbc];
      if (sum - min_sum == max_sum) {
        res.left_top_col = min_ltc;
        res.right_bot_col = rbc;
        break;
      }
      if (sum < min_sum) {
        min_sum = sum;
        min_ltc = rbc + 1;
      }   
    }
  }
  res.square = (res.right_bot_row - res.left_top_row + 1) *
               (res.right_bot_col - res.left_top_col + 1);
  *result = res;
}

static void * kadane_alg_threaded (void * pa)
{
  ARGS *pargs = (ARGS*) pa;
  cpu_set_t CPU_ID;
  int i;
  int *pps;

  CPU_ZERO(&CPU_ID);
  CPU_SET(pargs->CPU_ID, &CPU_ID);
  sched_setaffinity(0, sizeof(cpu_set_t), &CPU_ID);

  generate_matrix(pargs->n, pargs->m,
                  pargs->thread_number, pargs->total_threads, 
                  pargs->ps, pargs->len, pargs->prelen,
                  pargs->cycle, pargs->precycle, pargs->sync);  

  kadane_thread(pargs->n, pargs->m, pargs->MAX_ROW,
                pargs->total_threads, pargs->thread_number,
                pargs->ps, pargs->workspace, &(pargs->result));
                  
  return 0;
}

int kadane_alg (
  int n,                          // rows in matrix
  int m,                          // columns in matrix
  int total_threads,              // number of threads
  int start_thread,
  int active_threads,
  int* ps,                        // (n x m) sums matrix
  int* workspace,                 // (max(MIN_LEN, M) + 2 * M) workspace
  answer *result,                 // result
  int A, int B, int M, int SEED0  // pseudo random
  )
{
  const int MIN_SUM = -20001;
  pthread_t *threads;
  ARGS *args;
  pthread_barrier_t sync;         // synchronize barrier
  int i;
  int *cycle, *precycle, *checked;
  int len, prelen, min_len, MAX_ROW;
  int div, rem;
  cpu_set_t affinity;
  
  CPU_ZERO(&affinity);
  for (i = 0; i < total_threads; ++i) { 
    CPU_SET((start_thread + i) % active_threads, &affinity);
  }
  sched_setaffinity(0, sizeof(cpu_set_t), &affinity);
  
  checked = workspace;
  precycle = checked + M;
  cycle = precycle + M;

  pthread_barrier_init(&sync, 0, total_threads);
  while (!(args = (ARGS*) malloc(total_threads * sizeof(ARGS))))
  {
    printf("Not enough memory for args, waiting...\n");
  }
  
  while (!(threads = (pthread_t*) malloc(total_threads * sizeof(pthread_t))))
  {
    printf("Not enough memory for threads, waiting...\n");
  }
  
  min_len = generate_cycle(n, m, &prelen, &len, 
                           checked, cycle, precycle, A, B, M, SEED0);
                           
  MAX_ROW = min(n, M);
  memset(checked, 0, min_len * sizeof(int));
  div = (prelen + min_len + m - 1);
  div /= m;
  if (min_len != 0) {
    rem = (div * m - prelen) % min_len;
    while (checked[rem] == 0) {
      checked[rem] = 1;
      ++div;
      rem = (rem + m) % min_len;
    }
  } else {
    ++div;
  }
  MAX_ROW = min(MAX_ROW, div);
                           
  for (i = 0; i < total_threads; ++i) {
    args[i].n = n;
    args[i].m = m;
    args[i].MAX_ROW = MAX_ROW;
    args[i].total_threads = total_threads;
    args[i].thread_number = i;
    args[i].CPU_ID = (start_thread + i) % active_threads;
    while (!(args[i].workspace = 
             (int*) malloc((2 * BLOCK_ROW * BLOCK_COL + BLOCK_ROW +
                            BLOCK_ROW * BLOCK_ROW) * sizeof(int)))) 
    {
      printf("Not enough memory for workspace, waiting...\n");                                   
    }
    
    args[i].len = len;
    args[i].prelen = prelen;
    args[i].cycle = cycle;
    args[i].precycle = precycle;
    args[i].sync = &sync;
    args[i].ps = ps;
  }
  
  for (i = 0; i < total_threads; ++i) {
    while (pthread_create (threads + i, NULL, 
                           kadane_alg_threaded,
                           args + i)) 
    {
      printf("Cannot create thread: %d\n, waiting...", i);
    }
  }
  
  result->sum = MIN_SUM;  
  for (i = 0; i < total_threads; ++i) {
    while (pthread_join(threads[i], NULL)) {
      printf("Cannot join thread: %d...\n", i);
    }
    free(args[i].workspace);
    if (args[i].result.sum > result->sum) {
      *result = args[i].result;
    }
  }
  
  free(ps);
  pthread_barrier_destroy(&sync);
  free(args);
  free(threads); 
  
  return 0;
}
