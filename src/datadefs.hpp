#ifndef DATADEFS_HPP
#define DATADEFS_HPP

#include<cstdlib>
#include<vector>
#include<set>
#include<string>
#include<math.h>
#include<map>
#include<cassert>
#include<iostream>
#include<algorithm>

using namespace std;

namespace datadefs
{
  
  //Numerical data type
  typedef float num_t;

  //NaN as represented by the program
  extern const num_t NUM_NAN;
  extern const num_t EPS;
  extern const num_t A;
  extern const num_t PI;
  extern const num_t LOG_OF_MAX_NUM;

  //NaNs supported in the delimited file 
  typedef string NAN_t;
  extern const set<NAN_t> NANs;

  void strv2catv(vector<string>& strvec, vector<num_t>& catvec);
  void strv2numv(vector<string>& strvec, vector<num_t>& numvec);

  num_t str2num(string& str);

  inline bool isNAN(const string& str)
  {
    set<string>::const_iterator it(NANs.find(str));
    if(it == NANs.end()) { return(false); } else { return(true); } 
  }

  inline bool isNAN(const num_t value)
  {
    if(value != value) //== datadefs::NUM_NAN) 
      { return(true); } else { return(false); }
  }

 
  inline bool isNAN(const vector<num_t>& data)
  {
    for(size_t i = 0; i < data.size(); ++i)
      {
	if( data[i] != data[i]) //datadefs::isNAN(data[i]))
	  { return(true); }
      }
    return(false);
  }

  
  //DEPRECATED
  //void findNANs(vector<num_t>& data, vector<size_t>& NANIcs);

  void meanVals(vector<num_t> const& data, num_t& mu, size_t& nreal);

  // Newly added function set which is to be used in calculating leaf node predictions
  //void mean(const vector<num_t>& data, num_t& mu, const size_t numClasses);
  //void mode(const vector<num_t>& data, num_t& mode, const size_t numClasses);
  //void gamma(const vector<num_t>& data, num_t& gamma, const size_t numClasses);
  void cardinality(const vector<num_t>& data, size_t& cardinality);

  void sqerr(vector<num_t> const& data, 
	     num_t& mu, 
	     num_t& se,
	     size_t& nreal);      

  void countRealValues(vector<num_t> const& data, size_t& nRealValues);

  //DEPRECATED
  void zerotrim(vector<num_t>& data);

  
  inline void forward_sqerr(const num_t x_n,
			    size_t& n,
			    num_t& mu,
			    num_t& se)
  {  
    if(x_n != x_n) { return; } 
    ++n;
    num_t mu_old = mu;
    mu += (x_n - mu) / n;
    //If there are already at least two data points, squared error can be calculated, otherwise assign se_left := 0.0
    if(n > 1) { se += (x_n - mu) * (x_n - mu_old);} else { se = 0.0; }
  } 
    
  inline void forward_backward_sqerr(const num_t x_n,
			      size_t& n_left,
			      num_t& mu_left,
			      num_t& se_left,
			      size_t& n_right,
			      num_t& mu_right,
			      num_t& se_right)
  {
    if( x_n != x_n ) { return; }
    assert(n_right > 0); ++n_left; --n_right;
    //As long as there are at least two data points on the "right" branch, squared error can be calculated, otherwise assign se_right := 0.0
    if(n_right > 1) {
	  datadefs::num_t mu_old = mu_right;
	  mu_right -= (x_n - mu_right) / n_right;
	  se_right -= (x_n - mu_right) * (x_n - mu_old); }
      else if(n_right == 1) {
	  mu_right -= (x_n - mu_right) / n_right;
	  se_right = 0.0; }
      else {
	  mu_right = 0.0;
	  se_right = 0.0; }
      //Add x_n to "left" and update mean and squared error
      datadefs::num_t mu_old = mu_left;
      mu_left += (x_n - mu_left) / n_left;
      //If there are already at least two data points on the "left" branch, squared error can be calculated, otherwise assign se_left := 0.0
      if(n_left > 1) { se_left += (x_n - mu_left) * (x_n - mu_old); } else { se_left = 0.0; }
  }
  
  void count_freq(vector<num_t> const& data, map<num_t,size_t>& cat2freq, size_t& nRealValues);

  void map_data(vector<num_t>& data, 
		map<num_t,vector<size_t> >& datamap,
		size_t& nRealValues);

  void gini(vector<num_t>& data,
	    num_t& giniIndex,
	    size_t& nRealValues);

  void gini(map<num_t,size_t>& cat2freq,
            num_t& giniIndex);

  void sqfreq(vector<num_t> const& data, 
	      map<num_t,size_t>& freq, 
	      size_t& sqFreq, 
	      size_t& nRealValues);

  void forward_sqfreq(const num_t x_n,
		      size_t& n,
		      map<num_t,size_t>& freq,
		      size_t& sqFreq);
  
  void forward_backward_sqfreq(const num_t x_n,
			       size_t& n_left,
			       map<num_t,size_t>& freq_left,
			       size_t& sf_left,
			       size_t& n_right,
			       map<num_t,size_t>& freq_right,
			       size_t& sf_right);
  
  void range(vector<size_t>& ics);

  //DEPRECATED
  void ttest(vector<num_t> const& x, 
	     vector<num_t> const& y, 
	     num_t& pvalue);

  void utest(vector<num_t> const& x,
	     vector<num_t> const& y,
	     num_t& pvalue);

  num_t erf(num_t x);

  //DEPRECATED
  void regularized_betainc(const num_t x, 
			   const size_t a, 
			   num_t& ibval);

  void spearman_correlation(vector<num_t> const& x, 
			    vector<num_t> const& y,
			    num_t& corr);

  void pearson_correlation(vector<num_t> const& x,
                            vector<num_t> const& y,
                            num_t& corr);


  //DEPRECATED
  void percentile(vector<num_t> x, const num_t alpha, num_t& prc);

  
  //A comparator functor that can be passed to STL::sort. Assumes that one is comparing first elements of pairs, first type being num_t and second T
  template <typename T> struct increasingOrder {
    bool operator ()(pair<datadefs::num_t,T> const& a, pair<datadefs::num_t,T> const& b)
    {
      return(a.first < b.first);
    }
  }; 

  //A comparator functor that can be passed to STL::sort. Assumes that one is comparing first elements of pairs, first type being num_t and second T
  template <typename T> struct decreasingOrder {
    bool operator ()(pair<datadefs::num_t,T> const& a, pair<datadefs::num_t,T> const& b)
    {
      return(a.first > b.first);
    }
  };
  
  struct freqIncreasingOrder {
    bool operator ()(const map<datadefs::num_t,size_t>::value_type& a, const map<datadefs::num_t,size_t>::value_type& b)
    {
      return(a.second < b.second);
    }
  };
  
  template <typename T1,typename T2> void make_pairedv(vector<T1> const& v1,
                                                       vector<T2> const& v2,
                                                       vector<pair<T1,T2> >& p)
  {
    assert(v1.size() == v2.size() && v2.size() == p.size());
    for(size_t i = 0; i < p.size(); ++i)
      {
	p[i] = make_pair(v1[i],v2[i]);
      }
  }


  template <typename T1,typename T2> void separate_pairedv(vector<pair<T1,T2> > const& p,
                                                           vector<T1>& v1,
                                                           vector<T2>& v2)
  {
    assert(v1.size() == v2.size() && v2.size() == p.size());
    for(size_t i = 0; i < p.size(); ++i)
      {
	v1[i] = p[i].first;
	v2[i] = p[i].second;
      }
  }

  void sortDataAndMakeRef(const bool isIncreasingOrder, vector<num_t>& data, vector<size_t>& refIcs);

  //Sorts a given input data vector of type T based on a given reference ordering of type vector<int>
  template <typename T> void sortFromRef(vector<T>& data, vector<size_t> const& refIcs)
  {
    assert(data.size() == refIcs.size());  
    vector<T> foo = data;
    int n = data.size();
    for (int i = 0; i < n; ++i)
      {
	data[i] = foo[refIcs[i]];
      }
  }
}

#endif
