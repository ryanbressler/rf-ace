#ifndef RFACE_HPP
#define RFACE_HPP

#include <cppunit/extensions/HelperMacros.h>
//#include "datadefs.hpp"
//#include "hash.hpp"
#include "rf_ace.hpp"

class RFACETest : public CppUnit::TestFixture {
  CPPUNIT_TEST_SUITE( RFACETest );
  CPPUNIT_TEST( test_filter );
  CPPUNIT_TEST( test_trainRF );
  CPPUNIT_TEST( test_trainGBT );;
  CPPUNIT_TEST( test_test );
  CPPUNIT_TEST_SUITE_END();

public:
  void setUp();
  void tearDown();
  
  void test_filter();
  void test_trainRF();
  void test_trainGBT();
  void test_test();

private:

  Treedata* filterData_;
  Treedata* trainData_;
  Treedata* testData_;
  
  RFACE* rface_;

};

void RFACETest::setUp() {

  bool useContrasts = true;
  filterData_ = new Treedata("test_103by300_mixed_nan_matrix.afm",'\t',':',useContrasts);
  trainData_ = new Treedata("test_103by300_mixed_nan_matrix.afm",'\t',':');
  testData_ = new Treedata("test_103by300_mixed_nan_matrix.afm",'\t',':');

  rface_ = new RFACE;

}
void RFACETest::tearDown() {

  delete filterData_;
  delete trainData_;
  delete testData_;

  delete rface_;

}

void RFACETest::test_filter() {

  size_t targetIdx = 0;

  vector<num_t> featureWeights(filterData_->nFeatures(),1.0);
  featureWeights[0] = 0;

  ForestOptions forestOptions;

  forestOptions.nTrees = 10;
  forestOptions.mTry = 10;
  forestOptions.nMaxLeaves = 100;

  FilterOptions filterOptions;

  RFACE::FilterOutput filterOutput = rface_->filter(filterData_,targetIdx,featureWeights,&forestOptions,&filterOptions);

}

void RFACETest::test_trainRF() {
  
  vector<num_t> featureWeights(trainData_->nFeatures(),1.0);
  featureWeights[0] = 0;

  ForestOptions forestOptions;

  forestOptions.forestType = ForestOptions::ForestType::RF;
  forestOptions.nTrees = 10;
  forestOptions.mTry = 10;
  forestOptions.nMaxLeaves = 100;

  for ( size_t targetIdx = 0; targetIdx < trainData_->nFeatures(); ++targetIdx ) {
    vector<num_t> featureWeights(trainData_->nFeatures(),1.0);
    featureWeights[targetIdx] = 0;
    rface_->train(trainData_,targetIdx,featureWeights,&forestOptions);
  }

}

void RFACETest::test_trainGBT() {

  ForestOptions forestOptions;

  forestOptions.forestType = ForestOptions::ForestType::GBT;
  forestOptions.nTrees = 10;
  forestOptions.mTry = 10;
  forestOptions.nMaxLeaves = 100;

  for ( size_t targetIdx = 0; targetIdx < trainData_->nFeatures(); ++targetIdx ) {
    vector<num_t> featureWeights(trainData_->nFeatures(),1.0);
    featureWeights[targetIdx] = 0;
    rface_->train(trainData_,targetIdx,featureWeights,&forestOptions);
  }

}

void RFACETest::test_test() {

}


// Registers the fixture into the test 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION( RFACETest ); 

#endif // RFACE_HPP
