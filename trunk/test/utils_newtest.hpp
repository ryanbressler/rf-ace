#ifndef UTILS_NEWTEST_HPP
#define UTILS_NEWTEST_HPP

#include <cstdlib>
#include <vector>
#include <map>
#include <cmath>

#include "datadefs.hpp"
#include "newtest.hpp"

using namespace std;
using datadefs::num_t;


void utils_newtest_categoricalFeatureSplitsNumericalTarget();
void utils_newtest_categoricalFeatureSplitsCategoricalTarget();

void utils_newtest() {

  newtest("Testing categorical feature splitter with numerical targets", &utils_newtest_categoricalFeatureSplitsNumericalTarget);
  newtest("Testing categorical feature splitter with categorical targets", &utils_newtest_categoricalFeatureSplitsCategoricalTarget);

}

void utils_newtest_categoricalFeatureSplitsNumericalTarget() {

  vector<num_t> fv = {1,1,1,2,2,2,3,3,3,4,4,4};
  vector<num_t> tv = {1,1,1,2,3,4,5,6,7,8,9,10};

  map<num_t,vector<size_t> > fmap_left,fmap_right;

  num_t DI = utils::categoricalFeatureSplitsNumericalTarget2(tv,fv,1,{1,2,3,4},fmap_left,fmap_right);
  
  num_t DI_ref = math::deltaImpurity_regr(math::mean(tv),12,math::mean({1,1,1,2,3,4}),6,math::mean({5,6,7,8,9,10}),6);
  
  newassert( fabs( DI - DI_ref ) < datadefs::EPS );

  fv = {1,1,1,1,1,1,1,1,1,1,1,1};

  DI = utils::categoricalFeatureSplitsNumericalTarget2(tv,fv,1,{1},fmap_left,fmap_right);

  DI_ref = 0;

  newassert( fabs( DI - DI_ref ) < datadefs::EPS );  
  
}

void utils_newtest_categoricalFeatureSplitsCategoricalTarget() {
  
  

}

#endif
