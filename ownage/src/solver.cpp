#include <string>
#include <algorithm>
#include <map>
#include <iostream>
#include <fstream>

#include "solver.h"
#include "prefix_tree.h"

using namespace::std;

typedef std::pair<int, int> suffix_pair;

#define MAX_STEPS 40
#define MAX_PART_LEN 9

static void parse_file (
  const char* filename,
  string_system* my_system,
  bool referal
  )
{
  ifstream stream (filename,ios::in);
  string name = "";
  string dna_sequence = "";

  while (!stream.eof ())
    {
      string new_str;
      getline (stream, new_str);
      if (new_str[0] == '>')
        {
          if (!name.empty ())
            {
              if (!referal)
                my_system->add_string (name, dna_sequence);
              else
                my_system->set_referal (dna_sequence);
              dna_sequence.clear ();
            }

          name = new_str.substr (1);
        }
      else
        {
          dna_sequence += new_str;
        }
    }
  if (!name.empty ())
    {
      if (!referal)
        my_system->add_string (name, dna_sequence);
      else
        my_system->set_referal (dna_sequence);
    }
}

void solve (
  size_t thread_index,
  size_t thread_count,
  int argc,
  char* argv[],
  string_system* my_system,
  pthread_barrier_t *sync
  )
{
  int min_match_length = atol (argv[2]);
  int part_len = min (min_match_length, MAX_PART_LEN);
  int shift = 1;
  for (int i = 1; i < part_len; ++i)
    shift *= 4;
  shift -= 1;

  for (int i = 3; i < argc; ++i)
    {
      if (i % thread_count == thread_index)
        {
          if (i == 3)
            parse_file (argv[i], my_system, true);
          else
            parse_file (argv[i], my_system, false);
        }
    }
  pthread_barrier_wait (sync);
  /// Reading done

  size_t length = my_system->ref.size ();
  const string &ref = my_system->ref;

  static size_t total_len = 0;

  if (thread_index == 0)
    {
      /// TODO: Make parallel
      sort (my_system->data.begin (), my_system->data.end ());
      for (size_t i = 0; i < my_system->data.size (); ++i)
        total_len += my_system->data[i].str.size ();
    }
  pthread_barrier_wait (sync);

  static prefix_tree *ref_tree = 0;
  static prefix_tree *str_tree = 0;

  if (length < total_len && thread_index == 0)
    {
      ref_tree = new prefix_tree (part_len);
      ref_tree->parse (ref);
    }

  pthread_barrier_wait (sync);

  if (ref_tree)
    {
      ref_tree->calculate_depth (thread_index, thread_count, sync);
    }

  pthread_barrier_wait (sync);
  /// Referal string now in memory
  size_t size = my_system->data.size ();

  for (size_t index = thread_index; index < size; index += thread_count)
    {
      const string &str = my_system->data[index].str;
      size_t str_length = str.size ();

      if (length > 100 && str_length > 100 && length * str_length > 1000000)
        continue;
      /// Applying this alg only on tiny strings
      size_t i, j;
      char new_char;
      matches_pair new_match;
      vector<matches_pair> matches;
      int* old_layer = new int [str_length + 1];
      int* new_layer = new int [str_length + 1];
      new_char = ref[0];
      for (j = 0; j < str_length; ++j)
        old_layer[j] = (str[j] == new_char);
      for (i = 1; i < length; ++i)
        {
          new_char = ref[i];
          new_layer[0] = (str[0] == new_char);
          for (j = 1; j < str_length; ++j)
            new_layer[j] = (old_layer[j - 1] + 1) * (str[j] == new_char);

          for (j = 0; j < str_length; ++j)
            {
              if (old_layer[j] >= min_match_length)
                {
                  //this match can be shifted on i and j
                  if (j + 1 < str_length && old_layer[j] <= new_layer[j + 1])
                    continue;
                  //this match can be shifted on j
                  if (j + 1 < str_length && old_layer[j] <= old_layer[j + 1])
                    continue;
                  //this match can be shifted on i
                  if (old_layer[j] <= new_layer[j])
                    continue;
                  new_match.ref_start = i - old_layer[j] + 1;
                  new_match.ref_end = i;
                  new_match.str_start = j - old_layer[j] + 2;
                  new_match.str_end = j + 1;
                  matches.push_back (new_match);
                }
            }
          swap (old_layer, new_layer);
        }
      i = length;
      for (j = 0; j < str_length; ++j)
        {
          if (old_layer[j] >= min_match_length)
            {
              //this match can be shifted on i and j
              /// FALSE
              //this match can be shifted on j
              if (j + 1 < str_length && old_layer[j] <= old_layer[j + 1])
                continue;
              //this match can be shifted on i
              /// FALSE
              new_match.ref_start = i - old_layer[j] + 1;
              new_match.ref_end = i;
              new_match.str_start = j - old_layer[j] + 2;
              new_match.str_end = j + 1;
              matches.push_back (new_match);
            }
        }
      FREE_ARRAY (old_layer);
      FREE_ARRAY (new_layer);
      my_system->data[index].str.clear ();
      my_system->add_matches (index, matches);
    }
  pthread_barrier_wait (sync);

  for (size_t index = 0; index < size; ++index)
    {
      const string &str = my_system->data[index].str;
      if (str.empty ())
        continue;
      size_t str_length = str.size ();


      map<int, int> *old_layer = new map<int, int> ();
      map<int, int> *new_layer = new map<int, int> ();
      map<int, int>::iterator it;

      vector<matches_pair> matches;
      matches_pair new_match;
      if (ref_tree)
        {
          int first = thread_index * str_length / thread_count;
          int last = (thread_index + 1) * str_length / thread_count;
          int moves_count = min (min_match_length, MAX_STEPS);
          int start_pos = max (first, moves_count - 1);
          int part_str = 0;
          for (int i = max (0, start_pos - part_len); i < start_pos; ++i)
            part_str = ((part_str & shift)<<2) + char_code (str[i]) - 1;

          int node, pos, len;
          bool can_continue = true;
          for (int i = start_pos; i <= last; ++i)
            {
              if ((size_t) i < str_length)
                {
                  part_str = ((part_str & shift)<<2) + char_code (str[i]) - 1;
                  if (ref_tree->part_index (part_str) < 0)
                    {
                      if (can_continue)
                        continue;
                    }
                  else
                    {
                      node = ref_tree->part_index (part_str);
                      pos  = ref_tree->part_pos (part_str);
                      for (len = part_len; node >= 0 && len < moves_count; ++len)
                        ref_tree->move (str[i - len], node, pos);

                      if (node >= 0 && len == moves_count)
                        {
                          vector<int> indexes;
                          ref_tree->get_indexes (node, indexes);
                          for (size_t k = 0; k < indexes.size (); ++k)
                            {
                              int j = indexes[k];
                              int result = len;
                              if (i == start_pos)
                                {
                                  while (i - result >= 0 && j - result >= 0 &&
                                         str[i - result] == ref[j - result])
                                    result++;
                                }
                              if (old_layer->count (j - 1))
                                result = old_layer->at (j - 1) + 1;
                              new_layer->insert (suffix_pair (j, result));
                            }
                        }
                    }
                }
              for (it = old_layer->begin (); it != old_layer->end (); ++it)
                {
                  int j = it->first;
                  int result = it->second;
                  if (result < min_match_length)
                    continue;
                  if (new_layer->count (j + 1))
                    continue;
                  if (old_layer->count (j + 1) && result <= old_layer->at (j + 1))
                    continue;
                  if (new_layer->count (j) && result <= new_layer->at (j))
                    continue;
                  new_match.ref_start = j - result + 2;
                  new_match.ref_end = j + 1;
                  new_match.str_start = i - result + 1;
                  new_match.str_end = i;
                  matches.push_back (new_match);
                }
              old_layer->clear ();
              swap (old_layer, new_layer);
              can_continue = old_layer->empty ();
            }
        }
      else
        {
          if (thread_index == 0)
            {
              str_tree = new prefix_tree (part_len);
              str_tree->parse (str);
            }
          pthread_barrier_wait (sync);

          str_tree->calculate_depth (thread_index, thread_count, sync);

          pthread_barrier_wait (sync);

          int first = thread_index * length / thread_count;
          int last = (thread_index + 1) * length / thread_count;
          int moves_count = min (min_match_length, MAX_STEPS);
          int start_pos = max (first, moves_count - 1);
          int part_str = 0;
          for (int i = max (0, start_pos - part_len); i < start_pos; ++i)
            part_str = ((part_str & shift)<<2) + char_code (ref[i]) - 1;

          int node, pos, len;
          bool can_continue = true;
          for (int i = start_pos; i <= last; ++i)
            {
              if ((size_t) i < length)
                {
                  part_str = ((part_str & shift)<<2) + char_code (ref[i]) - 1;
                  if (str_tree->part_index (part_str) < 0)
                    {
                      if (can_continue)
                        continue;
                    }
                  else
                    {
                      node =str_tree->part_index (part_str);
                      pos  = str_tree->part_pos (part_str);
                      for (len = part_len; node >= 0 && len < moves_count; ++len)
                        str_tree->move (ref[i - len], node, pos);

                      if (node >= 0 && len == moves_count)
                        {
                          vector<int> indexes;
                          str_tree->get_indexes (node, indexes);
                          for (size_t k = 0; k < indexes.size (); ++k)
                            {
                              int j = indexes[k];
                              int result = len;
                              if (i == start_pos)
                                {
                                  while (i - result >= 0 && j - result >= 0 &&
                                         ref[i - result] == str[j - result])
                                    result++;
                                }
                              if (old_layer->count (j - 1))
                                result = old_layer->at (j - 1) + 1;
                              new_layer->insert (suffix_pair (j, result));
                            }
                        }
                    }
                }
              for (it = old_layer->begin (); it != old_layer->end (); ++it)
                {
                  int j = it->first;
                  int result = it->second;
                  if (result < min_match_length)
                    continue;
                  if (new_layer->count (j + 1))
                    continue;
                  if (old_layer->count (j + 1) && result <= old_layer->at (j + 1))
                    continue;
                  if (new_layer->count (j) && result <= new_layer->at (j))
                    continue;
                  new_match.ref_start = i - result + 1;
                  new_match.ref_end = i;
                  new_match.str_start = j - result + 2;
                  new_match.str_end = j + 1;
                  matches.push_back (new_match);
                }
              old_layer->clear ();
              swap (old_layer, new_layer);
              can_continue = old_layer->empty ();
            }
        }
      sort (matches.begin (), matches.end ());
      my_system->add_matches (index, matches);
      pthread_barrier_wait (sync);
      /// Calculation done
      if (thread_index + 1 == thread_count)
        /// TODO: Make parallel, merge sort needed to combine matches
        sort (my_system->data[index].matches.begin (), my_system->data[index].matches.end ());

      if (thread_index == 0)
        {
          FREE (str_tree);
          FREE (old_layer);
          FREE (new_layer);
        }
    }
  if (thread_index == 0)
    FREE (ref_tree);
}
