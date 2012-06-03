#ifndef ARRAY_OP_H
#define ARRAY_OP_H

#define BLOCK_ROW 64
#define BLOCK_COL 64
#define MIN_LEN 256
#define MAX_GEN 4

typedef struct _answer {
  int left_top_row;
  int left_top_col;
  int right_bot_row; 
  int right_bot_col;
  int sum;
  int square;
} answer;

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
  );
  
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
  );
  
#endif /* ARRAY_OP_H */
