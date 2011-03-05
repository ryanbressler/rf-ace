#include "datadefs.hpp"
#include<math.h>
#include<cassert>
#include<sstream>
#include<algorithm>
#include<iostream>

using namespace std;

//const datadefs::cat_t datadefs::cat_nan = -1;
const datadefs::num_t datadefs::num_nan = sqrt(-1.0);

const string initNANs[] = {"NA","NAN"};
const set<datadefs::NAN_t> datadefs::NANs(initNANs,initNANs+2);

void toupper(string& str)
{
  int (*pf)(int) = toupper;
  transform(str.begin(), str.end(), str.begin(), pf);
}

void datadefs::strv2catv(vector<string>& strvec, vector<datadefs::num_t>& catvec, map<string,size_t>& str2valmap)
{
  assert(strvec.size() == catvec.size());

  //Reset the map
  map<string,size_t> foo;
  str2valmap = foo;
  size_t val(0);

  //Map unique strings to values and store values in catvec as floats 
  for(size_t i = 0; i < strvec.size(); ++i)
    {
      //Transform string to uppercase
      toupper(strvec[i]);

      //If the string is not defined to NaN
      if(!datadefs::is_nan(strvec[i]))
	{
	  map<string,size_t>::iterator it;

	  //Try to find the string in the map. If it's not found, add the map...
	  it = str2valmap.find(strvec[i]);
	  if(it == str2valmap.end())
	    {
	      str2valmap.insert(pair<string,size_t>(strvec[i],val));
	      ++val;
	    }
	  //...and use the map to set the value for the output vector
	  catvec[i] = float(str2valmap[strvec[i]]);
	}
      //If the string is defined to NaN, however...
      else
	{
	  catvec[i] = datadefs::num_nan;
	}
    }  
}

void datadefs::strv2numv(vector<string>& strvec, vector<datadefs::num_t>& numvec)
{
  assert(strvec.size() == numvec.size());
  
  for(size_t i = 0; i < numvec.size(); ++i)
    {
      toupper(strvec[i]);
      if(!datadefs::is_nan(strvec[i]))
	{
	  numvec[i] = str2num(strvec[i]);
	}
      else
	{
	  numvec[i] = datadefs::num_nan;
	}
    }
}

datadefs::num_t datadefs::str2num(string& str)
{
  stringstream ss(str);
  datadefs::num_t ret;
  ss >> ret;
  return(ret);
}

//datadefs::num_t datadefs::cat2num(datadefs::cat_t value)
//{
//  return(float(value));
//}

bool datadefs::is_nan(string& str)
{
  set<string>::iterator it(NANs.find(str));
  if(it == NANs.end())
    {
      return(false);
    }
  else
    {
      return(true);
    }
}
