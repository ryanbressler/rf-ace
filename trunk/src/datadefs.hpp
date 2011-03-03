#ifndef DATADEFS_HPP
#define DATADEFS_HPP

#include<cstdlib>
#include<vector>
#include<set>
#include<string>

using namespace std;

namespace datadefs
{
  
  typedef int cat_t;
  typedef float num_t;

  extern cat_t cat_nan;
  extern num_t num_nan;

  //typedef string NAN_t;
  extern set<string> NANs;
  //extern string NAN;

  num_t cat2num(cat_t value);

  bool is_nan(cat_t value);
  bool is_nan(num_t value);

}

#endif
