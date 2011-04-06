#include <cstdlib>
#include <cassert>
#include <string>
#include <iostream>
#include <map>
#include "datadefs.hpp"

using namespace std;
using datadefs::num_t;


int main()
{

  size_t n(6);
  vector<num_t> numdata(n);
  vector<num_t> catdata(n);
  
  size_t n_acc_num(0);
  size_t n_acc_cat(0);
  num_t mu_acc(0.0);
  num_t se_acc(0.0);
  map<num_t,size_t> freq_acc;
  num_t sf_acc(0.0);
  
  for(size_t i = 0; i < n; ++i)
    {
      num_t numval(1.0*i);
      numdata[i] = numval;
      datadefs::forward_sqerr(numval,n_acc_num,mu_acc,se_acc);
      cout << "forward:\tn=" << n_acc_num << "\tmu=" << mu_acc << "\tse=" << se_acc << endl;
    }
  num_t mu,se;
  size_t nreal;
  datadefs::sqerr(numdata,mu,se,nreal);

  cout << "=> total:\tn=" << nreal << "\tmu=" << mu << "\tse=" << se << endl;

  for(size_t i = 0; i < n; ++i)
    {
      num_t catval(i);
      catdata[i] = catval;
      datadefs::forward_sqfreq(catval,n_acc_cat,freq_acc,sf_acc);
      cout << "forward:\tn=" << n_acc_cat << "\tnfreq=" << freq_acc.size() << "\tsf=" << sf_acc << endl;
    }

  num_t sf;
  map<num_t,size_t> freq;
  datadefs::sqfreq(catdata,freq,sf,nreal);
  cout << " => total:\tn=" << nreal << "\tnfreq=" << freq.size() << "\tsf=" << sf << endl;

  

  
  return(EXIT_SUCCESS);
}
