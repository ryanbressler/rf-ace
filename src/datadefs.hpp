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
  
  //Categorical data type
  //typedef int cat_t;
  //Numerical data type
  typedef float num_t;

  //NaNs as represented by the program
  //extern const cat_t cat_nan;
  extern const num_t num_nan;

  //NaNs supported in the delimited file 
  typedef string NAN_t;
  extern const set<NAN_t> NANs;

  void strv2catv(vector<string>& strvec, vector<num_t>& catvec, map<string,size_t>& str2valmap);
  void strv2numv(vector<string>& strvec, vector<num_t>& numvec);

  //cat_t str2cat(string& str);
  num_t str2num(string& str);

  //Function to convert categorical data to numerical
  //num_t cat2num(cat_t value);

  bool is_nan(string& str);

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
    //bool operator ()(pair<datadefs::cat_t,T> const& a, pair<datadefs::num_t,T> const& b)
    //{
    //  if(a.first < b.first || b.first == -1)
    //    {
    //      return(true);
    //    }
    //  else
    //    {
    //      return(false);
    //    }
    //}
  }; 
}

#endif
