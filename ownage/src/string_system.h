#ifndef STRING_SYSTEM_H
#define STRING_SYSTEM_H
#include <pthread.h>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

using std::string;
using std::vector;

class matches_pair
{
 public:
  int ref_start, ref_end;
  int str_start, str_end;
  friend bool operator< (const matches_pair &a, const matches_pair &b)
  {
    return (a.ref_end < b.ref_end || (a.ref_end == b.ref_end && a.str_end < b.str_end));
  }
};

class string_data
{
 public:
  string_data (const string& name, const string& str) :
    name (name), str (str)
  {
  }
  string_data ()
  {
  }

  friend bool operator< (const string_data& a, const string_data& b)
  {
    return a.name < b.name;
  }

  vector<matches_pair> matches;
  string name, str;
};

class string_system
{
public:
  string_system ()
  {
  }
  void add_string (const string& name, const string& str);
  void print ();
  void set_referal (const string& str);
  void add_matches (const size_t index, const vector<matches_pair> &matches);

  string ref;
  vector<string_data> data;
};

#endif // STRING_SYSTEM_H
