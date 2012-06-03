#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <pthread.h>

#include "array_op.h"

static __inline__ int
max(int a, int b) {
  return (a > b ? a: b);
}

static __inline__ int
min(int a, int b) {
  return (a < b ? a: b);
}

static __inline__ int
PRNG (
  int seed, 
  int a, 
  int b, 
  int m
  )
{
  return (a * seed + b) % m;
}

// add array src to dst: dst = dst + src
static __inline__ void
add_array (size_t n,           // array size
           const int *src,     // input array
           int *dst            // output array
  )
{
  size_t i;
#pragma simd
#pragma vector always
  for (i = 0; i < n; ++i) 
  {
    dst[i] += src[i];
  }
}

/*
 * Generate cycle and mean value with a, b, m, seed0
 * returns mininal length of cycle
 */
int
generate_cycle(
  int N, int M,                       // size of matrix
  int *prelen_val,                    // length of precycle
  int *len_val,                       // length of cycle
  int *checked,                       // (m) array
  int *cycle,                         // max(MIN_LEN, m) array
  int *precycle,                      // (m) array
  int a, int b, int m, int seed0      // constants
  )
{
  int seed = seed0;
  int len, prelen;
  int mean, psum, remainder;
  int iter = 0;
  size_t dstSize = N * M;
  size_t i;
  long long sum;
  
  memset(checked, 0, m * sizeof(int));
  seed = PRNG(seed, a, b, m);
  while (checked[seed] == 0) {
    ++iter;
    checked[seed] = iter;
    seed = PRNG(seed, a, b, m);
  }
  prelen = checked[seed] - 1;
  len = iter - prelen;
  
  prelen = min(prelen, dstSize);  // no need to generate more than matrix
  seed = seed0;
  psum = 0;
  for (iter = 0; iter < prelen; ++iter) {
    seed = PRNG(seed, a, b, m);
    precycle[iter] = seed;
    if (iter < dstSize)
      psum += seed;
  }
  sum = (long long) psum;
  psum = 0;
  for (iter = 0; iter < len; ++iter) {
    seed = PRNG(seed, a, b, m);
    cycle[iter] = seed;
    psum += seed;
  }
  i = (dstSize - prelen);
  i /= len;
  sum += (long long) psum * i;
  i = prelen + i * len;
  for (; i < dstSize; ++i)
    sum += (long long) cycle[(i - prelen) % len];
  
  /* calculate the mean value. Avoid float logic when making rounding. */
  mean = (int) (sum / (long long) dstSize); /* updated line */
  remainder = (int) (sum % (long long) dstSize); /* updated line */
  mean += (remainder * 2 > (signed) dstSize) ? (1) : (0); /* updated line */
    
#pragma unroll(4)
  for (iter = 0; iter < prelen; ++iter)
    precycle[iter] -= mean;
#pragma unroll(4)
  for (iter = 0; iter < len; ++iter)
    cycle[iter] -= mean;
  for (iter = len; iter + len < MIN_LEN; iter += len) {
    memcpy(cycle + iter, cycle, len * sizeof(int));
  }
  *prelen_val = prelen;
  *len_val = iter;
  return len;
}

/*
 * Generate (N x M) matrix A using cycle and precycle
 */
void
generate_matrix (
  size_t N, size_t M,             // rows, columns
  size_t thread_number,
  size_t total_threads,
  int *matrix,                    // (n x m) matrix
  int len, int prelen,            // len and prelen of cycle
  int *cycle,                     // (len) cycle array
  int *precycle,                  // (prelen) array
  pthread_barrier_t *sync         // synchronize barrier
  )
{
  size_t first, last, i, j;
  size_t cycle_first;
  long long pos;
  int active_threads;
  int *pm;  // pointer on matrix
  
  active_threads = total_threads;
  total_threads = min(active_threads, MAX_GEN);
  if (thread_number < total_threads) {
  
    pos = (long long) N * M * thread_number;
    pos /= (long long) total_threads;
    first = pos;
    pos = (long long) N * M * (thread_number + 1);
    pos /= (long long) total_threads;
    last = pos;
    
    if (first < prelen) {
      memcpy(matrix + first, precycle + first, 
             (min(prelen, last) - first) * sizeof(int));
      first = min(prelen, last);  
    }
    if (first >= prelen) {
      while (first < last) {
        cycle_first = (first - prelen) % len;
        memcpy(matrix + first, cycle + cycle_first, 
               min(len - cycle_first, last - first) * sizeof(int));
        first +=min(len - cycle_first, last - first);
      }
    }
  
  }
  
  //////////////////////////////////////////////////////////////////////////////
  pthread_barrier_wait(sync);
  //////////////////////////////////////////////////////////////////////////////
  // Each thread must have at least 64 bytes or 16 ints
  total_threads = min(active_threads, (M + 15) / 16);
  
  if (thread_number < total_threads) {
    first = M * thread_number;
    first /= total_threads;
    last = M * (thread_number + 1);
    last /= total_threads;
    
    for (i = 1; i < N; ++i) {
      add_array(last - first, 
                matrix + (i - 1) * M + first,
                matrix + i * M + first);
    }
  }
  
  //////////////////////////////////////////////////////////////////////////////
  pthread_barrier_wait(sync);
  //////////////////////////////////////////////////////////////////////////////
  total_threads = active_threads;
  
  if (thread_number < total_threads) {
    first = N * thread_number;
    first /= total_threads;
    last = N * (thread_number + 1);
    last /= total_threads;
    
    for (i = first, pm = matrix + first * M; i < last; ++i, pm += M) {
      for (j = 1; j < M; ++j)
        pm[j] += pm[j - 1];
    }
  }
  
  //////////////////////////////////////////////////////////////////////////////
  pthread_barrier_wait(sync);
  //////////////////////////////////////////////////////////////////////////////
}
