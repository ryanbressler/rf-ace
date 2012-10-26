#include "distributions.hpp"

using datadefs::num_t;

distributions::RandInt::RandInt():
  rand_(0,datadefs::MAX_IDX) {
  this->seed( distributions::generateSeed() );
}

distributions::RandInt::RandInt(size_t seed):
  rand_(0,datadefs::MAX_IDX) {

  this->seed(seed);

}

distributions::RandInt::~RandInt() {
  
}

void distributions::RandInt::seed(size_t seed) {
  eng_.seed(seed);
}

num_t distributions::RandInt::uniform() {

  return( 1.0 * rand_(eng_) / datadefs::MAX_IDX );

}

distributions::PMF::PMF(const vector<num_t>& weights) {

  num_t sum = 0.0;

  for ( size_t i = 0; i < weights.size(); ++i ) {
    assert( weights[i] >= 0.0 );
    sum += weights[i];
  }

  num_t cumProb = 0.0;

  assert( sum > 0.0 );

  for ( size_t i = 0; i < weights.size(); ++i ) {
    if ( weights[i] > 0.0 ) {
      cumProb += weights[i] / sum;
      icdf_[cumProb] = i;
    }
  }

}

distributions::PMF::~PMF() {
  
}

size_t distributions::PMF::icdf(const num_t prob) {
  assert(prob >= 0.0 && prob < 1.0);

  return( icdf_.upper_bound(prob)->second );
}

void distributions::walkersalias(
		  int n,   // number of classes
		  double *p, // relative weights of each class
		  int nans,  // sample size to return
		  int *ans  // sample as an array of class indices 
                  ){

  int *a = (int*) calloc(n,sizeof(int));
  double *q, rU;
  int i,j,k;
  int *HL,*H,*L;
  HL = (int*) calloc(n,sizeof(int));
  q = (double*) calloc(n,sizeof(double));
  H = HL - 1; L = HL + n;
  double sum = 0;
  for(i = 0; i < n; i++)
    {
      sum += p[i];
    }
  for(i = 0; i < n; i++)
    {
      p[i] /= sum;
    }
  for(i = 0; i < n; i++)
    {
      q[i] = p[i] * n;
      if(q[i] < 1.) *++H = i; else *--L = i;
    }
  if(H >= HL && L < HL +n)
    {
      for(k = 0; k < n-1; k++)
	{
	  i = HL[k];
	  j = *L;
	  a[i] = j;
	  q[j] += q[i] - 1;
	  if(q[j] < 1.) L++;
	  if(L >= HL + n) break;
	}
    }
  for(i = 0; i < n; i++) q[i] += i;
  for(i = 0; i < nans; i++)
    {
      rU = (double) rand() / RAND_MAX * n;
      k = (int) rU;
      ans[i] = (rU < q[k]) ? k : a[k];
    }
  free(HL);
  free(q);
  free(a);
}
