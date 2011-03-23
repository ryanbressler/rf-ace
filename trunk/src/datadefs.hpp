#ifndef DATADEFS_HPP
#define DATADEFS_HPP

#include<cstdlib>
#include<vector>
#include<set>
#include<string>
#include<math.h>
#include<map>

using namespace std;

namespace datadefs
{
  
  //Numerical data type
  typedef float num_t;
  typedef num_t* num_tp;

  //NaN as represented by the program
  extern const num_t num_nan;

  //NaNs supported in the delimited file 
  typedef string NAN_t;
  extern const set<NAN_t> NANs;

  void strv2catv(vector<string>& strvec, vector<num_t>& catvec, map<string,size_t>& str2valmap);
  void strv2numv(vector<string>& strvec, vector<num_t>& numvec);

  num_t str2num(string& str);

  bool is_nan(string& str);
  bool is_nan(num_t value);
  
  void sqerr(vector<num_t> const& data, 
	     num_t& mu, 
	     num_t& se,
	     size_t& nreal);      
  
  void update_sqerr(const num_t x_n,
		    const size_t n_left,
		    num_t& mu_left,
		    num_t& se_left,
		    const size_t n_right,
		    num_t& mu_right,
		    num_t& se_right);

  void count_freq(vector<num_t>& data, map<num_t,size_t>& cat2freq);

  void map_data(vector<num_t>& data, 
		map<num_t,vector<size_t> >& datamap);

  void gini(vector<num_t>& data,
	    num_t& gi);

  void gini(map<num_t,size_t>& cat2freq,
            num_t& gi);

  //void gini(map<num_t,vector<size_t> >& datamap, 
  //	    num_t& gi);
  
  void update_gini(num_t x_n,
		   const size_t n_left,
		   map<num_t,size_t>& cat2freq_left,
		   num_t& gi_left,
		   const size_t n_right,
		   map<num_t,size_t>& cat2freq_right,
		   num_t& gi_right);

  //TODO: ALL SORTING TO BE PERFORMED WITH VECTORS OF POINTERS TO DATA (THIS WILL REDUCE THE AMOUNT OF COPYING)

  //A comparator functor that can be passed to STL::sort. Assumes that one is comparing first elements of pairs, first type being num_t and second T
  template <typename T> struct ordering {
    bool operator ()(pair<datadefs::num_t,T> const& a, pair<datadefs::num_t,T> const& b)
    {
      if(a.first < b.first || b.first != b.first)
	{
	  return(true);
	}
      else
	{
	  return(false);
	}      
    }
  }; 

  struct paired_ptr_ordering {
    bool operator ()(pair<datadefs::num_tp,datadefs::num_tp> const& a, pair<datadefs::num_tp,datadefs::num_tp> const& b)
    {
      if(*(a.first) < *(b.first) || *(b.first) != *(b.first))
        {
          return(true);
        }
      else
        {
          return(false);
        }
    }
  };


}

#endif
