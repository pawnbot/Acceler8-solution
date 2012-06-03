#include <iostream>
#include <cstdio>

#include "string_system.h"

using std::cout;

static pthread_mutex_t operation_mutex = PTHREAD_MUTEX_INITIALIZER;

void string_system::print ()
{
  for (size_t i = 0; i < data.size (); ++i)
    {
      cout << data[i].name << "\n";
      for (size_t j = 0; j < data[i].matches.size (); ++j)
        printf ("%d %d %d %d\n", data[i].matches[j].ref_start, data[i].matches[j].ref_end,
                data[i].matches[j].str_start, data[i].matches[j].str_end);
    }
}

void string_system::add_string (const string& name, const string& str)
{
  pthread_mutex_lock (&operation_mutex);
  data.push_back (string_data (name, str));
  pthread_mutex_unlock (&operation_mutex);
}

void string_system::set_referal (const string& str)
{
  pthread_mutex_lock (&operation_mutex);
  ref = str;
  pthread_mutex_unlock (&operation_mutex);
}

void string_system::add_matches (const size_t index, const vector<matches_pair> &matches)
{
  pthread_mutex_lock (&operation_mutex);
  data[index].matches.insert (data[index].matches.end (), matches.begin (), matches.end ());
  pthread_mutex_unlock (&operation_mutex);
}
