#ifndef DATADEFS_NEWTEST_HPP
#define DATADEFS_NEWTEST_HPP

#include <cstdlib>
#include <limits>

#include "newtest.hpp"
#include "datadefs.hpp"

void datadefs_newtest_isnan();
void datadefs_newtest_countRealValues();
void datadefs_newtest_mapData();
void datadefs_newtest_increasingOrder();
void datadefs_newtest_decreasingOrder();
void datadefs_newtest_freqIncreasingOrder();
void datadefs_newtest_makePairedV();
void datadefs_newtest_separatePairedV();

void datadefs_newtest() {

  newtest( "Testing interpretation of missing values", &datadefs_newtest_isnan );
  newtest( "Testing counting of real values", &datadefs_newtest_countRealValues );
  newtest( "Testing data mapping", &datadefs_newtest_mapData );
  newtest( "Testing incr. ordering operator", &datadefs_newtest_increasingOrder );
  newtest( "Testing decr. ordering operator", &datadefs_newtest_decreasingOrder );
  newtest( "Testing counting frequency", &datadefs_newtest_freqIncreasingOrder );
  newtest( "Testing pairing of a vector", &datadefs_newtest_makePairedV );
  newtest( "Testing separation of a pried vector", &datadefs_newtest_separatePairedV );

}

void datadefs_newtest_isnan() {
 
  newassert(datadefs::isNAN_STR("na"));
  newassert(datadefs::isNAN_STR("NA"));
  newassert(datadefs::isNAN_STR("nan"));
  newassert(datadefs::isNAN_STR("NaN"));
  newassert(datadefs::isNAN_STR("NAN"));
  newassert(datadefs::isNAN_STR("?"));
  
  newassert(!datadefs::isNAN_STR("2"));
  newassert(!datadefs::isNAN_STR("@data"));
  newassert(!datadefs::isNAN_STR("NAte"));
  newassert(!datadefs::isNAN_STR("name"));
  newassert(!datadefs::isNAN_STR("na?"));

  newassert(numeric_limits<datadefs::num_t>::has_quiet_NaN);
  newassert(datadefs::isNAN(datadefs::NUM_NAN));
  
  if (numeric_limits<datadefs::num_t>::has_infinity) {
    newassert(!datadefs::isNAN(numeric_limits<datadefs::num_t>::infinity()));
    newassert(!datadefs::isNAN(-numeric_limits<datadefs::num_t>::infinity()));
  }
  
  newassert(!datadefs::isNAN((datadefs::num_t)0.0));
  newassert(!datadefs::isNAN((datadefs::num_t)-1.0));
  newassert(!datadefs::isNAN((datadefs::num_t)1.0));
  
}

void datadefs_newtest_countRealValues() {
  
  vector<datadefs::num_t> data;
  size_t nRealValues = static_cast<size_t>(-1);
  for (int i = 0; i < 50; ++i) {
    data.push_back(static_cast<datadefs::num_t>(i));
  }
  
  datadefs::countRealValues(data, nRealValues);
  newassert(nRealValues == 50);
  
  // Our data vector is defined as const in our signature. Since we're not
  //  testing edge cases of non-trivial memory corruption, we ignore it here.
  
  // Interleave the original input with NaNs; verify we get the same results
  for (int i = 0; i < 50; ++i) {
    data.insert(data.begin() + (i*2), datadefs::NUM_NAN);
  }
  
  nRealValues = static_cast<size_t>(-1);
  
  datadefs::countRealValues(data, nRealValues);
  newassert(nRealValues == 50);
  
  // Ensure a vector containing only NaNs is handled as expected
  data.clear();
  for (int i = 0; i < 50; ++i) {
    data.push_back(datadefs::NUM_NAN);
  }
  
  nRealValues = static_cast<size_t>(-1);
  
  datadefs::countRealValues(data, nRealValues);
  newassert(nRealValues == 0);
}

void datadefs_newtest_mapData() {
  
  vector<datadefs::num_t> data;
  unordered_map<datadefs::num_t,vector<size_t> > datamap;
  size_t nRealValues;
  for (int i = 0; i < 50; ++i) {
    data.push_back(static_cast<datadefs::num_t>(i));
  }
  data.push_back(0.0);
  
  datadefs::map_data(data, datamap, nRealValues);
  newassert(datamap.size() == 50);
  newassert(nRealValues == 51);
  
  // Our data vector is defined as const in our signature. Since we're not
  //  testing edge cases of non-trivial memory corruption, we ignore it here.
  
  // Interleave the original input with NaNs; verify we get the same results
  for (int i = 0; i < 50; ++i) {
    data.insert(data.begin() + (i*2), datadefs::NUM_NAN);
  }
  
  nRealValues = static_cast<size_t>(-1);
  datadefs::map_data(data, datamap, nRealValues);
  newassert(datamap.size() == 50);
  newassert(nRealValues == 51);
  
  // Ensure a vector containing only NaNs is handled as expected
  data.clear();
  for (int i = 0; i < 50; ++i) {
    data.push_back(datadefs::NUM_NAN);
  }
  
  nRealValues = static_cast<size_t>(-1);
  datadefs::map_data(data, datamap, nRealValues);
  newassert(datamap.size() == 0);
  newassert(nRealValues == 0);
}

void datadefs_newtest_increasingOrder() {
  
  pair<datadefs::num_t, int> a(1.0, 2);
  pair<datadefs::num_t, int> b(2.0, 1);
  
  datadefs::increasingOrder<int> incOrder;
  newassert(incOrder.operator()(a,b));
  newassert(!incOrder.operator()(b,a));
  
  vector<pair<datadefs::num_t, int> > pairedTestVector;
  pairedTestVector.push_back(pair<datadefs::num_t, int>(5.0,6));
  pairedTestVector.push_back(pair<datadefs::num_t, int>(4.0,5));
  pairedTestVector.push_back(pair<datadefs::num_t, int>(3.0,4));
  pairedTestVector.push_back(pair<datadefs::num_t, int>(2.0,3));
  pairedTestVector.push_back(pair<datadefs::num_t, int>(1.0,2));
  pairedTestVector.push_back(pair<datadefs::num_t, int>(0.0,1));
  pairedTestVector.push_back(pair<datadefs::num_t, int>(-1.0,10));
  pairedTestVector.push_back(pair<datadefs::num_t, int>(
							datadefs::NUM_NAN,0));
  
  sort(pairedTestVector.begin(), pairedTestVector.end(), incOrder);
  newassert(pairedTestVector[0].first == -1.0);
  newassert(pairedTestVector[0].second == 10);
  for (int i = 1; i < 7; ++i) {
    newassert(pairedTestVector[i].first == static_cast<double>(i-1));
    newassert(pairedTestVector[i].second == i);
  }
  newassert(datadefs::isNAN(pairedTestVector[7].first));
  newassert(pairedTestVector[7].second == 0);
}

void datadefs_newtest_decreasingOrder() {
  
  pair<datadefs::num_t, int> a(1.0, 2);
  pair<datadefs::num_t, int> b(2.0, 1);
  
  datadefs::decreasingOrder<int> decOrder;
  newassert(decOrder.operator()(b,a));
  newassert(!decOrder.operator()(a,b));
  
  vector<pair<datadefs::num_t, int> > pairedTestVector;
  pairedTestVector.push_back(pair<datadefs::num_t, int>(5.0,6));
  pairedTestVector.push_back(pair<datadefs::num_t, int>(4.0,5));
  pairedTestVector.push_back(pair<datadefs::num_t, int>(3.0,4));
  pairedTestVector.push_back(pair<datadefs::num_t, int>(2.0,3));
  pairedTestVector.push_back(pair<datadefs::num_t, int>(1.0,2));
  pairedTestVector.push_back(pair<datadefs::num_t, int>(0.0,1));
  pairedTestVector.push_back(pair<datadefs::num_t, int>(-1.0,10));
  pairedTestVector.push_back(pair<datadefs::num_t, int>(
							datadefs::NUM_NAN,0));
  
  sort(pairedTestVector.begin(), pairedTestVector.end(), decOrder);
  for (int i = 0; i > 6; --i) {
    newassert(pairedTestVector[i].first == static_cast<double>(5-i));
    newassert(pairedTestVector[i].second == 5-i);
  }
  newassert(pairedTestVector[6].first == -1.0);
  newassert(pairedTestVector[6].second == 10);
  newassert(datadefs::isNAN(pairedTestVector[7].first));
  newassert(pairedTestVector[7].second == 0);
}

void datadefs_newtest_freqIncreasingOrder() {
  pair<datadefs::num_t, size_t> a(1.0, 2);
  pair<datadefs::num_t, size_t> b(2.0, 1);
  
  datadefs::freqIncreasingOrder freqIncOrder;
  newassert(freqIncOrder.operator()(b,a));
  newassert(!freqIncOrder.operator()(a,b));
  
  vector<pair<datadefs::num_t, size_t> > pairedTestVector;
  pairedTestVector.push_back(pair<datadefs::num_t, size_t>(5.0,6));
  pairedTestVector.push_back(pair<datadefs::num_t, size_t>(4.0,5));
  pairedTestVector.push_back(pair<datadefs::num_t, size_t>(3.0,4));
  pairedTestVector.push_back(pair<datadefs::num_t, size_t>(2.0,3));
  pairedTestVector.push_back(pair<datadefs::num_t, size_t>(1.0,2));
  pairedTestVector.push_back(pair<datadefs::num_t, size_t>(0.0,1));
  pairedTestVector.push_back(pair<datadefs::num_t, size_t>(-1.0,10));
  pairedTestVector.push_back(pair<datadefs::num_t, size_t>(datadefs::NUM_NAN,0));

  sort(pairedTestVector.begin(), pairedTestVector.end(), freqIncOrder);
  newassert(datadefs::isNAN(pairedTestVector[0].first));
  newassert(pairedTestVector[0].second == 0);
  for ( int i = 1; i < 7; ++i) {
    newassert(pairedTestVector[i].first == static_cast<double>(i-1));
    newassert(pairedTestVector[i].second == static_cast<size_t>(i));
  }
  newassert(pairedTestVector[7].first == -1.0);
  newassert(pairedTestVector[7].second == 10);
}

void datadefs_newtest_makePairedV() {
  vector<int> v1(50,1);
  vector<string> v2(50,"a");
  vector<pair<int,string> > p(50,pair<int,string>(0,""));
  
  datadefs::make_pairedv(v1,v2,p);
  for (int i = 0; i < 50; ++i) {
    newassert(v1[i] == 1);
    newassert(strcmp(v2[i].c_str(),"a") == 0);
    newassert(p[i].first == 1);
    newassert(strcmp(p[i].second.c_str(),"a") == 0);
  }
}

void datadefs_newtest_separatePairedV() {
  vector<pair<int,string> > p(50,pair<int,string>(1,"a"));
  vector<int> v1(50,0);
  vector<string> v2(50,"");

  datadefs::separate_pairedv(p,v1,v2);
  for (int i = 0; i < 50; ++i) {
    newassert(v1[i] == 1);
    newassert(strcmp(v2[i].c_str(),"a") == 0);
    newassert(p[i].first == 1);
    newassert(strcmp(p[i].second.c_str(),"a") == 0);
  }
}



#endif
