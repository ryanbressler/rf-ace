#include "datadefs.hpp"
#include<math.h>

using namespace std;

datadefs::cat_t datadefs::cat_nan = -1;
datadefs::num_t datadefs::num_nan = sqrt(-1.0);

bool datadefs::is_nan(datadefs::cat_t value)
{
  if(value == datadefs::cat_nan)
    {
      return(true);
    }
  return(false);
}

bool datadefs::is_nan(datadefs::num_t value)
{
  return(value != datadefs::num_nan);
}
