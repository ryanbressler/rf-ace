#include "datadefs.hpp"
#include<math.h>
#include<cassert>
#include<sstream>
#include<algorithm>
#include<iostream>

using namespace std;

const datadefs::cat_t datadefs::cat_nan = -1;
const datadefs::num_t datadefs::num_nan = sqrt(-1.0);

const string initNANs[] = {"NA","NAN"};
const set<datadefs::NAN_t> datadefs::NANs(initNANs,initNANs+2);

void toupper(string& str)
{
  int (*pf)(int) = toupper;
  transform(str.begin(), str.end(), str.begin(), pf);
}

void datadefs::strv2catv(vector<string>& strvec, vector<datadefs::cat_t>& catvec)
{
  assert(strvec.size() == catvec.size());
  
  for(size_t i = 0; i < strvec.size(); ++i)
    {
      toupper(strvec[i]);
      catvec[i] = str2cat(strvec[i]);
    }  
}

void datadefs::strv2numv(vector<string>& strvec, vector<datadefs::num_t>& numvec)
{
  assert(strvec.size() == numvec.size());
  
  for(size_t i = 0; i < numvec.size(); ++i)
    {
      toupper(strvec[i]);
      numvec[i] = str2num(strvec[i]);
    }
}

datadefs::cat_t datadefs::str2cat(string& str)
{
  stringstream ss(str);
  string foo;
  datadefs::cat_t ret;
  ss >> foo;
  assert(ss.eof());
  set<string>::iterator it(NANs.find(foo));
  if(it == NANs.end())
    {
      ss.clear();
      ss.str("");
      ss << foo;
      ss >> ret;
      assert(ret > datadefs::cat_nan);
      return(ret);
    }
  else
    {
      return(datadefs::cat_nan);
    }
}

datadefs::num_t datadefs::str2num(string& str)
{
  stringstream ss(str);
  string foo;
  datadefs::num_t ret;
  ss >> foo;
  assert(ss.eof());
  set<string>::iterator it(NANs.find(foo));
    if(it == NANs.end())
    {
      ss.clear();
      ss.str("");
      ss << foo;
      ss >> ret;
      return(ret);
    }
  else
    {
      return(datadefs::num_nan);
    }
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
