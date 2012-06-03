#ifndef PREFIX_TREE_H
#define PREFIX_TREE_H
#include <pthread.h>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

using std::vector;
using std::string;

const int ALPHABET_SIZE = 5;
const int INFINITY = 2000000000;

#define FREE(A) do { if (A) {delete (A);} (A) = 0; } while (0)
#define FREE_ARRAY(A) do { if (A) {delete[] (A);} (A) = 0; } while (0)

static inline int
char_code (char symbol)
{
  switch (symbol)
    {
      case ('A'):
        return 1;
      case ('C'):
        return 2;
      case ('G'):
        return 3;
      case ('T'):
        return 4;
      default:
        return symbol;
    }
}

/// prefix inverse tree
class prefix_tree
{
private:
  int max_nodes;
  int total_len;
  int len;
  int *str;
  int **transitions;
  int *start, *end, *parent, *buf;
  int node_id, node_pos, node_link;
  int next_link;
  void append (int new_char);
  void cleanup ();
  int *str_index;
  int *str_pos;
  int part_len;
  int total_parts;

public:
  explicit prefix_tree (int part_len);
  ~prefix_tree ();
  void parse (const string& new_str);
  void calculate_depth (int thread_index, int thread_count, pthread_barrier_t *sync);
  void copy_indexes (int *part_index);
  inline int edge_len (int node_id)
  {
    return end[node_id] - start[node_id] + 1;
  }

  inline int part_index (int index)
  {
    return str_index[index];
  }

  inline int part_pos (int index)
  {
    return str_pos[index];
  }

  void inline get_indexes (int node, vector<int>& result)
  {
    vector<int> queue;
    queue.push_back (node);
    int v, i;
    for (size_t pos = 0; pos < queue.size (); ++pos)
      {
        v = queue[pos];
        if (end[v] == INFINITY && buf[v] > 0)
          {
            result.push_back (buf[v] - 1);
          }

        for (i = 0; i < ALPHABET_SIZE; ++i)
          if (transitions[v][i] != -1)
            {
              queue.push_back (transitions[v][i]);
            }
      }
  }

  inline void move (char symbol, int &node, int &pos)
  {
    int ch = char_code (symbol);
    if (end[node] < pos)
      {
        node = transitions[node][ch];
        if (node >= 0)
          pos = start[node] + 1;
        else
          return;
      }
    else
      {
        if (str[pos] == ch)
          pos++;
        else
          {
            node = -1;
          }
      }
  }
};

#endif // PREFIX_TREE_H
