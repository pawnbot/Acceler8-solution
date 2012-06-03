#include "prefix_tree.h"

prefix_tree::prefix_tree (int part_len) : part_len (part_len)
{
  transitions = 0;
  start = 0;
  end = 0;
  parent = 0;
  buf = 0;
  str = 0;
  str_index = 0;
  str_pos = 0;
  total_parts = 1;
  for (int i = 0; i < part_len; ++i)
    total_parts *= 4;
}

void prefix_tree::cleanup ()
{
  if (transitions)
    {
      for (int i = 0; i < max_nodes; ++i)
        {
          FREE_ARRAY (transitions[i]);
        }
      FREE_ARRAY (transitions);
    }
  FREE_ARRAY (start);
  FREE_ARRAY (end);
  FREE_ARRAY (parent);
  FREE_ARRAY (buf);
  FREE_ARRAY (str);
  FREE_ARRAY (str_index);
  FREE_ARRAY (str_pos);
}

prefix_tree::~prefix_tree ()
{
  cleanup ();
}

/// Ukkonen tree algorithm without recursion
inline void prefix_tree::append (int new_char)
{
  bool place_found = false;
  while (!place_found)
    {
      if (end[node_id] < node_pos)
        {
          if (transitions[node_id][new_char] >= 0)
            {
              node_id = transitions[node_id][new_char];
              node_pos = start[node_id];
            }
          else
            {
              transitions[node_id][new_char] = node_link;
              start[node_link] = len;
              parent[node_link] = node_id;
              node_id = buf[node_id];
              node_pos = end[node_id] + 1;
              node_link++;
              continue;
            }
        }
      if (node_pos >= 0 && new_char != str[node_pos])
        {
          next_link = node_link + 1;
          start[next_link] = len;
          parent[next_link] = node_link;
          start[node_link] = start[node_id];
          end[node_link] = node_pos - 1;
          parent[node_link] = parent[node_id];
          transitions[node_link][new_char] = next_link;
          transitions[node_link][str[node_pos]] = node_id;
          start[node_id] = node_pos;
          parent[node_id] = node_link;
          transitions[parent[node_link]][str[start[node_link]]] = node_link;
          node_id = buf[parent[node_link]];
          node_pos = start[node_link];
          while (node_pos <= end[node_link])
            {
              node_id   = transitions[node_id][str[node_pos]];
              node_pos += edge_len (node_id);
            }
          next_link++;
          if (node_pos == end[node_link] + 1)
            buf[node_link] = node_id;
          else
            buf[node_link] = next_link;
          node_pos = end[node_id] + 1 - (node_pos - end[node_link] - 1);
          node_link = next_link;
          continue;
        }
      else
        {
          node_pos++;
          place_found = true;
        }
    }
}

void prefix_tree::parse (const string &new_str)
{
  cleanup ();
  total_len = new_str.size () + 1;
  max_nodes = total_len * 2 + 1;
  transitions = new int* [max_nodes];
  start = new int [max_nodes];
  end = new int [max_nodes];
  parent = new int [max_nodes];
  /// buf array will store suffix links before Ukkonen and depth after
  buf = new int [max_nodes];
  for (int i = 0; i < max_nodes; ++i)
    {
      end[i] = INFINITY;
      transitions[i] = new int [ALPHABET_SIZE];
      for (int j = 0; j < ALPHABET_SIZE; ++j)
        transitions[i][j] = -1;
    }
  /// Adding Joker and Root
  buf[0] = 1;
  start[0] = -1;
  end[0] = -1;
  start[1] = -1;
  end[1] = -1;
  /// All links from Joker to Root
  for (int j = 0; j < ALPHABET_SIZE; ++j)
    transitions[1][j] = 0;
  /// Default position
  node_link = 2;
  node_id = 0;
  node_pos = 0;
  str = new int [total_len];
  for (int i = 0; i < (total_len - 1); ++i)
    {
      str[i] = char_code (new_str[(total_len - 1) - i - 1]);
    }
  str[total_len - 1] = (char) 0;
  for (len = 0; len < total_len; ++len)
    append (str[len]);

  FREE_ARRAY (parent);

  memset (buf, 0, max_nodes * sizeof (int));
  /// Calculating positions for all strings with length = part_len
  str_index = new int [total_parts];
  str_pos = new int [total_parts];
  memset (str_index, -1, total_parts * sizeof (int));
  vector<int> queue;
  vector<int> position;
  vector<int> str_part;
  queue.push_back (0);
  str_part.push_back (0);
  position.push_back (0);
  int new_v, new_pos, index, i;
  size_t new_len = 1;
  int length = 0;
  int shift = 1;
  for (size_t pos = 0; pos < queue.size (); ++pos)
    {
      if (pos == new_len)
        {
          new_len = queue.size ();
          length++;
          shift = (shift<<2);
        }
      if (length == part_len)
        break;

      for (i = 1, index = str_part[pos]; i < ALPHABET_SIZE; ++i, index += shift)
        {
          new_v = queue[pos];
          new_pos = position[pos];
          move (i, new_v, new_pos);
          if (new_v >= 0)
            {
              queue.push_back (new_v);
              position.push_back (new_pos);
              str_part.push_back (index);
              if (length + 1 == part_len)
                {
                  str_index[index] = new_v;
                  str_pos[index] = new_pos;
                }
            }
        }
    }
}

void prefix_tree::calculate_depth (int thread_index, int thread_count, pthread_barrier_t *sync)
{
  /// buf array stores depth of every node
  /// TODO: Make parallel
  const int solo_size = 1024;
  static int shared[solo_size + ALPHABET_SIZE + 1];
  static int length = 0;
  static int shared_pos = 0;
  int v, i;
  if (thread_index == 0)
    {
      shared[0] = 0;
      length = 1;
      for (shared_pos = 0; shared_pos < length && length < solo_size; ++shared_pos)
        {
          v = shared[shared_pos];
          if (end[v] == INFINITY)
            {
              buf[v] += (total_len - 2) - start[v];
            }
          else
            buf[v] += edge_len (v);

          for (i = 0; i < ALPHABET_SIZE; ++i)
            if (transitions[v][i] != -1)
              {
                shared[length] = transitions[v][i];
                buf[transitions[v][i]] = buf[v];
                length++;
              }
        }
    }
  pthread_barrier_wait (sync);
  vector<int> queue;
  for (int pos = shared_pos + thread_index; pos < length; pos += thread_count)
    queue.push_back (shared[pos]);
  for (size_t pos = 0; pos < queue.size (); ++pos)
    {
      v = queue[pos];
      if (end[v] == INFINITY)
        {
          buf[v] += (total_len - 2) - start[v];
        }
      else
        buf[v] += edge_len (v);

      for (i = 0; i < ALPHABET_SIZE; ++i)
        if (transitions[v][i] != -1)
          {
            queue.push_back (transitions[v][i]);
            buf[transitions[v][i]] = buf[v];
          }
    }
}

void prefix_tree::copy_indexes (int *part_index)
{
  memcpy (part_index, str_index, total_parts * sizeof (int));
}
