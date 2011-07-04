#include "datadefs.hpp"
#include<math.h>
#include<cassert>
#include<sstream>
#include<algorithm>
#include<iostream>
#include<limits>

using namespace std;

const datadefs::num_t datadefs::NUM_NAN = numeric_limits<float>::infinity();
const datadefs::num_t datadefs::EPS = 1e-12;
const datadefs::num_t datadefs::PI = 3.1415926535;
const datadefs::num_t datadefs::A = 0.140012;
const datadefs::num_t datadefs::LOG_OF_MAX_NUM = 70.0;

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
  size_t val = 0;

  //Map unique strings to values and store values in catvec as floats 
  for(size_t strIdx = 0; strIdx < strvec.size(); ++strIdx)
    {
      //Transform string to uppercase
      toupper(strvec[strIdx]);

      //If the string is not defined to NaN
      if(!datadefs::isNAN(strvec[strIdx]))
	{
	  map<string,size_t>::iterator it;

	  //Try to find the string in the map. If it's not found, add the map...
	  it = str2valmap.find(strvec[strIdx]);
	  if(it == str2valmap.end())
	    {
	      str2valmap.insert(pair<string,size_t>(strvec[strIdx],val));
	      ++val;
	    }
	  //...and use the map to set the value for the output vector
	  catvec[strIdx] = float(str2valmap[strvec[strIdx]]);
	}
      //If the string is defined to NaN, however...
      else
	{
	  catvec[strIdx] = datadefs::NUM_NAN;
	}
    }  
}

void datadefs::strv2numv(vector<string>& strvec, vector<datadefs::num_t>& numvec)
{
  assert(strvec.size() == numvec.size());
  
  for(size_t strIdx = 0; strIdx < strvec.size(); ++strIdx)
    {
      toupper(strvec[strIdx]);
      if(!datadefs::isNAN(strvec[strIdx]))
	{
	  numvec[strIdx] = str2num(strvec[strIdx]);
	}
      else
	{
	  numvec[strIdx] = datadefs::NUM_NAN;
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


bool datadefs::isNAN(const string& str)
{
  set<string>::const_iterator it(NANs.find(str));
  if(it == NANs.end())
    {
      return(false);
    }
  else
    {
      return(true);
    }
}



bool datadefs::isNAN(const datadefs::num_t value)
{
  if(value == datadefs::NUM_NAN)
    {
      return(true);
    }
  else
    {
      return(false);
    }
}

bool datadefs::isNAN(const vector<num_t>& data)
{
  for(size_t i = 0; i < data.size(); ++i)
    {
      if(datadefs::isNAN(data[i]))
	{
	  return(true);
	}
    }
  return(false);
}

//This is a very poorly designed function
void datadefs::findNANs(vector<num_t>& data, vector<size_t>& NANIcs)
{
  NANIcs.clear();
  size_t n = data.size();
  for(size_t i = 0; i < n; ++i)
    {
      if(datadefs::isNAN(data[i]))
        {
          //data.erase(data.begin() + i);
	  NANIcs.push_back(i);
        }
    }
}

void datadefs::mean(vector<datadefs::num_t> const& data, datadefs::num_t& mu, size_t& nRealValues)
{
 
  nRealValues = 0;
  mu = 0.0;
 
  for(size_t i = 0; i < data.size(); ++i)
    {
      if(!datadefs::isNAN(data[i]))
        {
          ++nRealValues;
          mu += data[i];
        }
    }

  if(nRealValues > 0)
    {
      mu /= nRealValues;
    }

}

void datadefs::sqerr(vector<datadefs::num_t> const& data, 
		     datadefs::num_t& mu, 
		     datadefs::num_t& se,
		     size_t& nRealValues)
{
    
  datadefs::mean(data,mu,nRealValues);
  
  se = 0.0;
  for(size_t i = 0; i < data.size(); ++i)
    {
      if(!datadefs::isNAN(data[i]))
	{
	  se += pow(data[i] - mu,2);
	}
    }
}

void datadefs::countRealValues(vector<num_t> const& data, size_t& nRealValues)
{
  nRealValues = 0;
  for(size_t i = 0; i < data.size(); ++i)
    {
      if(!datadefs::isNAN(data[i]))
	{
	  ++nRealValues;
	}
    }  
}

//DEPRECATED
void datadefs::zerotrim(vector<datadefs::num_t>& data)
{
  int n = data.size();
  for(int i = n-1; i >= 0; --i)
    {
      if(fabs(data[i]) < datadefs::EPS)
	{
	  data.erase(data.begin() + i);
	}
    }
}

void datadefs::sortDataAndMakeRef(vector<num_t>& data, vector<size_t>& refIcs)
{
  //cout << "sort_and_make_ref: in the beginning" << endl;
  //assert(v.size() == ref_ics.size());
  vector<pair<num_t,size_t> > pairedv(data.size());
  refIcs.resize(data.size());
  datadefs::range(refIcs);
  //cout << "sort_and_make_ref: used range()" << endl;
  datadefs::make_pairedv<num_t,size_t>(data,refIcs,pairedv);
  //cout << "sort_and_make_ref: made pairedv" << endl;
  sort(pairedv.begin(),pairedv.end(),datadefs::ordering<num_t>());
  //cout << "sort_and_make_ref: pairedv sorted" << endl;
  datadefs::separate_pairedv<num_t,size_t>(pairedv,data,refIcs);
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

  if(datadefs::isNAN(x_n))
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
		    datadefs::num_t& giniIndex,
		    size_t& nRealValues)
{
  map<datadefs::num_t,size_t> freq;
  datadefs::count_freq(data,freq,nRealValues);
  datadefs::gini(freq,giniIndex);
}

void datadefs::gini(map<datadefs::num_t,size_t>& cat2freq, 
		    datadefs::num_t& giniIndex)
{
  giniIndex = 0.0;
  size_t n = 0;
  map<datadefs::num_t,size_t>::const_iterator it;
  for(it = cat2freq.begin(); it != cat2freq.end(); ++it)
    {
      size_t freq_new = it->second;
      giniIndex += freq_new * freq_new;
      n += freq_new;
    }
  if(n)
    {
      giniIndex = 1 - giniIndex / ( n*n );
    }
}

void datadefs::count_freq(vector<datadefs::num_t> const& data, map<datadefs::num_t,size_t>& cat2freq, size_t& nRealValues)
{
  cat2freq.clear();
  map<datadefs::num_t,size_t>::const_iterator it;
  nRealValues = 0;
  for(size_t i = 0; i < data.size(); ++i)
    {
      if(!datadefs::isNAN(data[i]))
	{
	  ++nRealValues;
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
			size_t& nRealValues)
{
  datamap.clear();
  map<datadefs::num_t,vector<size_t> >::iterator it;
  nRealValues = 0;
  for(size_t i = 0; i < data.size(); ++i)
    {
      if(!datadefs::isNAN(data[i]))
	{
	  ++nRealValues;
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
		      size_t& sqFreq, 
		      size_t& nRealValues)
{
  sqFreq = 0;
  datadefs::count_freq(data,freq,nRealValues);
  for(map<datadefs::num_t,size_t>::const_iterator it(freq.begin()); it != freq.end(); ++it)
    {
      size_t freq_new = it->second;
      sqFreq += freq_new * freq_new;
    }
}

void datadefs::forward_sqfreq(const datadefs::num_t x_n,
			      size_t& n,
			      map<datadefs::num_t,size_t>& freq,
			      size_t& sqFreq)
{

  if(datadefs::isNAN(x_n))
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
      sqFreq += 1;
      //cout << "  sf_new=" << sf << endl;

      //If not, add a new category and set its frequency to 1...
      freq.insert(pair<datadefs::num_t,size_t>(x_n,1));

    }
  else
    {
      sqFreq += 2*freq[x_n] + 1;
      ++freq[x_n];
    }
  //cout << "  sf_new=" << sf << endl;

}


void datadefs::forward_backward_sqfreq(const datadefs::num_t x_n,
				       size_t& n_left,
				       map<datadefs::num_t,size_t>& freq_left, 
				       size_t& sf_left,
				       size_t& n_right,
				       map<datadefs::num_t,size_t>& freq_right,
				       size_t& sf_right)
{
 
  if(datadefs::isNAN(x_n))
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

//DEPRECATED
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

    if(var_x < datadefs::EPS && var_y < datadefs::EPS)
      {
	if(fabs(mean_x - mean_y) < datadefs::EPS)
	  {
	    pvalue = 0.5;
	  }
	else if(mean_x > mean_y)
	  {
	    pvalue = 0.0;
	  }
	else
	  {
	    pvalue = 1.0;
	  }
	return;
      }

    size_t v;
    datadefs::num_t sp,tvalue,ttrans;
    
    v = nreal_x + nreal_y - 2;
    sp = sqrt(((nreal_x-1) * var_x + (nreal_y-1) * var_y) / v);
    tvalue = (mean_x - mean_y) / (sp * sqrt(1.0 / nreal_x + 1.0 / nreal_y));
    
    ttrans = (tvalue+sqrt(pow(tvalue,2) + v)) / (2 * sqrt(pow(tvalue,2) + v));

    datadefs::regularized_betainc(ttrans,nreal_x - 1,pvalue);	
    pvalue = 1-pvalue;
    
}

//DEPRECATED
void datadefs::regularized_betainc(const datadefs::num_t x,
                                   const size_t a,
                                   datadefs::num_t& ibval)
{

  ibval = 0.0;
  
  datadefs::num_t jfac = 1;
  for(size_t i = 1; i < a; ++i)
    {
      jfac += log(static_cast<datadefs::num_t>(i));
    }
  
  datadefs::num_t kfac = 1;
  for(size_t i = a+1; i < 2*a; ++i)
    {
      kfac += log(static_cast<datadefs::num_t>(i));
    }
  
  for(size_t i = a; i < 2*a; ++i)
    {
      jfac += log(static_cast<datadefs::num_t>(i));
      kfac += log(static_cast<datadefs::num_t>(2*a - i));
      datadefs::num_t temp = kfac - jfac + i*log(x) + (2*a-1-i)*log(1-x);
      ibval += exp(temp);
      
      //cout << jfac << "\t" << kfac << "\t" << ibval << endl;
    }

  if(ibval > 1.0)
    {
      ibval = 1.0;
    }
  
}

void datadefs::utest(vector<datadefs::num_t> const& x,
		     vector<datadefs::num_t> const& y,
		     datadefs::num_t& pvalue)
{
  
  num_t uvalue = 0.0;
  size_t m = x.size();
  size_t n = y.size();

  for(size_t i = 0; i < m; ++i)
    {

      bool xnan = datadefs::isNAN(x[i]);

      for(size_t j = 0; j < n; ++j)
	{
	 
	  bool ynan = datadefs::isNAN(y[j]);

	  if(!xnan && ynan)
	    {
	      uvalue += 1;
	    }
	  else if(!xnan && !ynan && x[i] > y[j])
	    {
	      uvalue += 1;
	    }
  
	}
    }
  

  num_t mu = 1.0*m*n/2.0;
  num_t s = sqrt(1.0*m*n*(n+m+1)/12.0);

  //cout << uvalue << " " << mu << " " << s;

  pvalue = 1.0 - 0.5 * ( 1 + datadefs::erf( (uvalue-mu) / (s*sqrt(2.0)) ) );
  //cout << " ==> " << pvalue << endl;

}

datadefs::num_t datadefs::erf(datadefs::num_t x)
{  
  
  num_t x2 = x*x;

  num_t sgn;
  if(x < 0.0)
    {
      sgn = -1.0;
    }
  else
    {
      sgn = 1.0;
    }

  return( sgn*sqrt(1.0 - exp(-x2*(4.0/datadefs::PI+datadefs::A*x2) / (1+datadefs::A*x2))) ); 

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
      if(!datadefs::isNAN(x[i]) && !datadefs::isNAN(y[i]))
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



//DEPRECATED
void datadefs::percentile(vector<datadefs::num_t> x, 
			  const datadefs::num_t alpha, 
			  datadefs::num_t& prc)
{

  sort(x.begin(),x.end());
  
  datadefs::num_t k((x.size()-1) * alpha);
  num_t f = floor(k);
  num_t c = ceil(k);
	  
  if(fabs(f - c) < datadefs::EPS)
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


