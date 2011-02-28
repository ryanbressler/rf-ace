#ifndef DATADEFS_HPP
#define DATADEFS_HPP

#include<cstdlib>

using namespace std;

namespace datadefs
{
  
  typedef int cat_t;
  typedef float num_t;

  extern cat_t cat_nan;
  extern num_t num_nan;

  num_t cat2num(cat_t value);

  bool is_nan(cat_t value);
  bool is_nan(num_t value);

}

#endif
