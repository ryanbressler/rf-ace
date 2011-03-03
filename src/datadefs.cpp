#include "datadefs.hpp"
#include<math.h>

using namespace std;

const datadefs::cat_t datadefs::cat_nan = -1;
const datadefs::num_t datadefs::num_nan = sqrt(-1.0);

const string initNANs[] = {"NA","NAN"};
const set<datadefs::NAN_t> datadefs::NANs(initNANs,initNANs+2);

void datadefs::str2cat(vector<string>& strvec, vector<datadefs::cat_t>& catvec)
{
  
}

void datadefs::str2num(vector<string>& strvec, vector<datadefs::num_t>& numvec)
{

}

datadefs::num_t datadefs::cat2num(datadefs::cat_t value)
{
  return(float(value));
}

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
