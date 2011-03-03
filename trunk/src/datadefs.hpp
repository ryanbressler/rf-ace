#ifndef DATADEFS_HPP
#define DATADEFS_HPP

#include<cstdlib>
#include<vector>
#include<set>
#include<string>

using namespace std;

namespace datadefs
{
  
  //Categorical data type
  typedef int cat_t;
  //Numerical data type
  typedef float num_t;

  //NaNs as represented by the program
  extern const cat_t cat_nan;
  extern const num_t num_nan;

  //NaNs supported in the delimited file 
  typedef string NAN_t;
  extern const set<NAN_t> NANs;

  void strv2catv(vector<string>& strvec, vector<cat_t>& catvec);
  void strv2numv(vector<string>& strvec, vector<num_t>& numvec);

  cat_t str2cat(string& str);
  num_t str2num(string& str);

  //Function to convert categorical data to numerical
  num_t cat2num(cat_t value);

  bool is_nan(cat_t value);
  bool is_nan(num_t value);

}

#endif
