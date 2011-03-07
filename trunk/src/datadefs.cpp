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

bool datadefs::is_nan(datadefs::num_t value)
{
  if(value != value)
    {
      return(true);
    }
  else
    {
      return(false);
    }
}

void datadefs::sqerr(vector<datadefs::num_t> const& data, 
		     vector<size_t> const& ics,
		     size_t& n,
		     datadefs::num_t& mu, 
		     datadefs::num_t& se)
{
  
  
  se = 0.0;
  if(datadefs::is_nan(data[ics[0]]))
    {
      mu = 0.0;
      n = 0;
    }
  else
    {
      mu = data[ics[0]];
      n = 1;
    }
  size_t ntot = ics.size();
  
  for(size_t i = 1; i < ntot; ++i)
    {
      if(!datadefs::is_nan(data[ics[i]]))
	{
	  ++n;
	  mu += data[ics[i]];
	}
    }
  
  mu /= n;
  
  //This should be computed iteratively inside the previous loop (speed-up)!
  for(size_t i = 0; i < ntot; ++i)
    {
      if(!datadefs::is_nan(data[ics[i]]))
	{
	  se += pow(data[ics[i]] - mu,2);
	}
    }
}


//Assuming x_n is a current member of the "right" branch, subtract it from "right" and add it to "left", and update the branch data counts, means, and squared errors. NOTE: NaN checks not implemented
void datadefs::update_sqerr(const datadefs::num_t x_n,
			    const size_t n_left,
			    datadefs::num_t& mu_left,
			    datadefs::num_t& se_left,
			    const size_t n_right,
			    datadefs::num_t& mu_right,
			    datadefs::num_t& se_right,
			    datadefs::num_t& mu_old)
{
  
  //Subtract x_n from "right" and update mean and squared error
  mu_old = mu_right;
  mu_right -= (x_n - mu_right) / n_right;

  //As long as there are at least two data points on the "right" branch, squared error can be calculated, otherwise assign se_right := 0.0
  if(n_right > 1)
    {
      se_right -= (x_n - mu_right) * (x_n - mu_old);
    }
  else
    {
      se_right = 0.0;
    }

  //Add x_n to "left" and update mean and squared error
  mu_old = mu_left;
  mu_left += (x_n - mu_left) / n_left;

  //If there are already at least two data points on the "left" branch, squared error can be calculated, otherwise assign se_left := 0.0
  if(n_left > 1)
    {
      se_left += (x_n - mu_left) * (x_n - mu_old);
    }
  else
    {
      se_left = 0.0;
    }
}
