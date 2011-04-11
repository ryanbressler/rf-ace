#include "datadefs.hpp"
#include<math.h>
#include<cassert>
#include<sstream>
#include<algorithm>
#include<iostream>
#include<limits>

//#include <boost/math/distributions/students_t.hpp>

using namespace std;

//const datadefs::cat_t datadefs::cat_nan = -1;
const datadefs::num_t datadefs::num_nan = numeric_limits<float>::infinity();
const datadefs::num_t datadefs::eps = 1e-12;

const string initNANs[] = {"NA","NAN"};
const set<datadefs::NAN_t> datadefs::NANs(initNANs,initNANs+2);

void toupper(string& str)
{
  int (*pf)(int) = toupper;
  transform(str.begin(), str.end(), str.begin(), pf);
}

void datadefs::strv2catv(vector<string>& strvec, vector<datadefs::num_t>& catvec)
{
  assert(strvec.size() == catvec.size());

  map<string,size_t> str2valmap;

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

bool datadefs::is_nan(const string& str)
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

bool datadefs::is_nan(const datadefs::num_t value)
{
  if(value == datadefs::num_nan)
    {
      return(true);
    }
  else
    {
      return(false);
    }
}

void datadefs::sqerr(vector<datadefs::num_t> const& data, 
		     datadefs::num_t& mu, 
		     datadefs::num_t& se,
		     size_t& nreal)
{
  
  nreal = 0;
  mu = 0.0;
  se = 0.0;
  
  size_t n(data.size());

  if(n == 0)
    {
      return;
    }

  if(!datadefs::is_nan(data[0]))
    {
      mu = data[0];
      nreal = 1;
    }
  
  for(size_t i = 1; i < n; ++i)
    {
      if(!datadefs::is_nan(data[i]))
	{
	  ++nreal;
	  mu += data[i];
	}
    }
  
  if(nreal > 0)
    {
      mu /= nreal;
    }
  
  //This should be computed iteratively inside the previous loop (speed-up)!
  for(size_t i = 0; i < n; ++i)
    {
      if(!datadefs::is_nan(data[i]))
	{
	  se += pow(data[i] - mu,2);
	}
    }
}

/*
  void datadefs::sqerr2(vector<datadefs::num_t> const& x, 
  vector<datadefs::num_t> const& y,  
  datadefs::num_t& se,
  size_t& nreal)
  {
  
  mu_x = 0.0;
  mu_y = 0.0;
  se = 0.0;
  
  
  size_t n(data.size());
  
  if(n == 0)
  {
  return;
  }
  
  if(!datadefs::is_nan(data[0]))
  {
  mu = data[0];
  nreal = 1;
  }
  
  for(size_t i = 1; i < n; ++i)
  {
  if(!datadefs::is_nan(data[i]))
  {
  ++nreal;
  mu += data[i];
  }
  }
  
  if(nreal > 0)
  {
  mu /= nreal;
  }
  
  //This should be computed iteratively inside the previous loop (speed-up)!
  for(size_t i = 0; i < n; ++i)
  {
  if(!datadefs::is_nan(data[i]))
  {
  se += pow(data[i] - mu,2);
  }
  }
  }
*/

void datadefs::count_real_values(vector<num_t> const& data, size_t& nreal)
{
  nreal = 0;
  for(size_t i = 0; i < data.size(); ++i)
    {
      if(!datadefs::is_nan(data[i]))
	{
	  ++nreal;
	}
    }  
}

void datadefs::forward_sqerr(const datadefs::num_t x_n,
			     size_t& n,
			     datadefs::num_t& mu,
			     datadefs::num_t& se)
{
  if(datadefs::is_nan(x_n))
    {
      return;
    }

  ++n;
  
  datadefs::num_t mu_old(mu);

  mu += (x_n - mu) / n;

  //If there are already at least two data points, squared error can be calculated, otherwise assign se_left := 0.0
  if(n > 1)
    {
      se += (x_n - mu) * (x_n - mu_old);
    }
  else
    {
      se = 0.0;
    }


}

//Assuming x_n is a current member of the "right" branch, subtract it from "right" and add it to "left", and update the branch data counts, means, and squared errors. NOTE: NaN checks not implemented
void datadefs::forward_backward_sqerr(const datadefs::num_t x_n,
				      size_t& n_left,
				      datadefs::num_t& mu_left,
				      datadefs::num_t& se_left,
				      size_t& n_right,
				      datadefs::num_t& mu_right,
				      datadefs::num_t& se_right)
{

  if(datadefs::is_nan(x_n))
    {
      return;
    }
  
  assert(n_right > 0);

  ++n_left;
  --n_right;

  //As long as there are at least two data points on the "right" branch, squared error can be calculated, otherwise assign se_right := 0.0
  if(n_right > 1)
    {
      datadefs::num_t mu_old(mu_right);
      mu_right -= (x_n - mu_right) / n_right;
      se_right -= (x_n - mu_right) * (x_n - mu_old);
    }
  else if(n_right == 1)
    {
      mu_right -= (x_n - mu_right) / n_right;
      se_right = 0.0;
    }
  else
    {
      mu_right = 0.0;
      se_right = 0.0;
    }

  //Add x_n to "left" and update mean and squared error
  datadefs::num_t mu_old = mu_left;
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

void datadefs::gini(vector<datadefs::num_t>& data,
		    datadefs::num_t& gi,
		    size_t& nreal)
{
  map<datadefs::num_t,size_t> freq;
  datadefs::count_freq(data,freq,nreal);
  datadefs::gini(freq,gi);
}

void datadefs::gini(map<datadefs::num_t,size_t>& cat2freq, 
		    datadefs::num_t& gi)
{
  gi = 0.0;
  size_t n(0);
  map<datadefs::num_t,size_t>::const_iterator it;
  for(it = cat2freq.begin(); it != cat2freq.end(); ++it)
    {
      gi += pow(it->second,2);
      n += it->second;
    }
  if(n)
    {
      gi = 1-gi/pow(n,2);
    }
}

void datadefs::count_freq(vector<datadefs::num_t> const& data, map<datadefs::num_t,size_t>& cat2freq, size_t& nreal)
{
  cat2freq.clear();
  map<datadefs::num_t,size_t>::const_iterator it;
  nreal = 0;
  for(size_t i = 0; i < data.size(); ++i)
    {
      if(!datadefs::is_nan(data[i]))
	{
	  ++nreal;
	  it = cat2freq.find(data[i]);
	  if(it == cat2freq.end())
	    {
	      cat2freq.insert(pair<datadefs::num_t,size_t>(data[i],1));
	    }
	  else
	    {
	      ++cat2freq[data[i]];
	    }
	}
    }
}


void datadefs::map_data(vector<datadefs::num_t>& data, 
			map<datadefs::num_t,vector<size_t> >& datamap, 
			size_t& nreal)
{
  datamap.clear();
  map<datadefs::num_t,vector<size_t> >::iterator it;
  nreal = 0;
  for(size_t i = 0; i < data.size(); ++i)
    {
      if(!datadefs::is_nan(data[i]))
	{
	  ++nreal;
	  it = datamap.find(data[i]);
	  if(it == datamap.end())
	    {
	      vector<size_t> foo;
	      foo.push_back(i);
	      datamap.insert(pair<datadefs::num_t,vector<size_t> >(data[i],foo));
	    }
	  else
	    {
	      it->second.push_back(i);
	    }
	}
    }
}

void datadefs::sqfreq(vector<datadefs::num_t> const& data, 
		      map<datadefs::num_t,size_t>& freq, 
		      datadefs::num_t& sf, 
		      size_t& nreal)
{
  sf = 0.0;
  datadefs::count_freq(data,freq,nreal);
  for(map<datadefs::num_t,size_t>::const_iterator it(freq.begin()); it != freq.end(); ++it)
    {
      sf += pow(it->second,2);
    }
}

void datadefs::forward_sqfreq(const datadefs::num_t x_n,
			      size_t& n,
			      map<datadefs::num_t,size_t>& freq,
			      datadefs::num_t& sf)
{

  if(datadefs::is_nan(x_n))
    {
      return;
    }

  ++n;
  //cout << "sf_old=" << sf;
  //Check if the value already exists
  map<datadefs::num_t,size_t>::const_iterator it(freq.find(x_n));
  if(it == freq.end())
    {
      //cout << "sf_old=" << sf;
      sf += 1;
      //cout << "  sf_new=" << sf << endl;

      //If not, add a new category and set its frequency to 1...
      freq.insert(pair<datadefs::num_t,size_t>(x_n,1));

    }
  else
    {
      sf += 2*freq[x_n] + 1;
      ++freq[x_n];
    }
  //cout << "  sf_new=" << sf << endl;

}


void datadefs::forward_backward_sqfreq(const datadefs::num_t x_n,
				       size_t& n_left,
				       map<datadefs::num_t,size_t>& freq_left, 
				       datadefs::num_t& sf_left,
				       size_t& n_right,
				       map<datadefs::num_t,size_t>& freq_right,
				       datadefs::num_t& sf_right)
{
 
  if(datadefs::is_nan(x_n))
    {
      return;
    }

  assert(n_right > 0);

  ++n_left;
  --n_right;

  //Check if the value already exists on left
  map<datadefs::num_t,size_t>::const_iterator it(freq_left.find(x_n));
  if(it == freq_left.end())
    {
      sf_left += 1;

      //If not, add a new category and set its frequency to 1...
      freq_left.insert(pair<datadefs::num_t,size_t>(x_n,1));
      
    }
  else
    {
      sf_left += 2*freq_left[x_n] + 1;
      ++freq_left[x_n];
    }

  it = freq_right.find(x_n);
  assert(it != freq_right.end());
  assert(freq_right[x_n] > 0);

  sf_right -= 2*freq_right[x_n] - 1;
  --freq_right[x_n];

  if(freq_right[x_n] == 0)
    {
      freq_right.erase(x_n);
    }  
}


void datadefs::range(vector<size_t>& ics)
{
  for(size_t i = 0; i < ics.size(); ++i)
    {
      ics[i] = i;
    }
}

void datadefs::ttest(vector<datadefs::num_t> const& x, 
		     vector<datadefs::num_t> const& y, 
		     datadefs::num_t& pvalue)
{
    // Sample mean and variance of x
    datadefs::num_t mean_x = 0;
    datadefs::num_t var_x = 0;
    size_t nreal_x = 0;

    datadefs::sqerr(x, mean_x, var_x, nreal_x);

    assert(nreal_x > 1);

    var_x /= (nreal_x - 1);

    // Sample mean and variance of y
    datadefs::num_t mean_y = 0;
    datadefs::num_t var_y = 0;
    size_t nreal_y = 0;

    datadefs::sqerr(y, mean_y, var_y, nreal_y);
    
    assert(nreal_y > 1);

    var_y /= (nreal_y - 1);
    
    assert(nreal_x == nreal_y);

    if((fabs(mean_x - mean_y) < datadefs::eps && var_x < datadefs::eps && var_y < datadefs::eps) || var_x < datadefs::eps)
      {
        pvalue = 1;
	return;
      }

    size_t v;
    datadefs::num_t sp,t,ttrans;
    if(fabs(mean_y) < datadefs::eps && var_y < datadefs::eps) //Reduce to one-sample t-test
      {
	v = nreal_x - 1;
	sp = sqrt(var_x/nreal_x);
	t = sqrt(1.0*nreal_x)*mean_x / sqrt(var_x);
      }
    else //Two-sample t-test
      {
	v = nreal_x + nreal_y - 2;
	sp = sqrt(((nreal_x-1) * var_x + (nreal_y-1) * var_y) / v);
	t = (mean_x - mean_y) / (sp * sqrt(1.0 / nreal_x + 1.0 / nreal_y));
      }
    ttrans = (t+sqrt(pow(t,2) + v)) / (2 * sqrt(pow(t,2) + v));

    datadefs::regularized_betainc(ttrans,nreal_x - 1,pvalue);	
    pvalue = 1-pvalue;
    
    cout << mean_x << "\t" << var_x << "\t" << mean_y << "\t" << var_y << "\t" << v << "\t" << sp << "\t" << t << "\t" << ttrans << "\t" << pvalue << endl;

}

void datadefs::regularized_betainc(const num_t x,
				   const size_t a,
				   num_t& ibval)
{

  ibval = 0;

  num_t jfac = 1;
  for(size_t i = 1; i < a; ++i)
    {
      jfac *= i;
    }
  
  num_t kfac = 1;
  for(size_t i = a+1; i < 2*a; ++i)
    {
      kfac *= i;
    }

  for(size_t i = a; i < 2*a; ++i)
    {
      jfac *= i;
      kfac *= 2*a - i;
      ibval += kfac/jfac*powf(x,i)*powf(1-x,2*a-1-i);

      //cout << jfac << "\t" << kfac << "\t" << ibval << endl;
    }

}

void datadefs::spearman_correlation(vector<datadefs::num_t> const& x,
				    vector<datadefs::num_t> const& y,
				    datadefs::num_t& corr)
{
  assert(false);
}

void datadefs::pearson_correlation(vector<datadefs::num_t> const& x,
				   vector<datadefs::num_t> const& y,
				   datadefs::num_t& corr)
{

  corr = 0.0;

  datadefs::num_t mu_x,se_x,mu_y,se_y;
  size_t nreal_x,nreal_y;

  vector<datadefs::num_t> x_real;
  vector<datadefs::num_t> y_real;

  size_t n = x.size();
  assert(n == y.size());
  
  for(size_t i = 0; i < n; ++i)
    {
      if(!datadefs::is_nan(x[i]) && !datadefs::is_nan(y[i]))
	{
	  x_real.push_back(x[i]);
	  y_real.push_back(y[i]);
	}
    }

  datadefs::sqerr(x_real,mu_x,se_x,nreal_x);
  datadefs::sqerr(y_real,mu_y,se_y,nreal_y);
  assert(nreal_x == nreal_y);
  
  for(size_t i = 0; i < nreal_x; ++i)
    {
      corr += ( x_real[i] - mu_x ) * ( y_real[i] - mu_y ); 
    }

  corr /= sqrt(se_x*se_y);

}



void datadefs::percentile(vector<datadefs::num_t> x, 
			  const datadefs::num_t alpha, 
			  datadefs::num_t& prc)
{

  sort(x.begin(),x.end());
  
  datadefs::num_t k((x.size()-1) * alpha);
  num_t f = floor(k);
  num_t c = ceil(k);
	  
  if(fabs(f - c) < datadefs::eps)
    {
      prc = x[static_cast<size_t>(k)];
    }
  else
    {
      num_t d0 = x[static_cast<size_t>(f)] * (c - k);
      num_t d1 = x[static_cast<size_t>(c)] * (k - f);
      prc = d0+d1;
    }
}

/*
  void datadefs::beta_symmetric(size_t n, datadefs::num_t& b)
  {
  b = gamma(2*n)/(gamma(n)*gamma(n));
  }
*/
