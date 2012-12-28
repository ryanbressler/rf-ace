#ifndef STOCHASTICFORESTTEST_HPP
#define STOCHASTICFORESTTEST_HPP

#include <cppunit/extensions/HelperMacros.h>
#include "datadefs.hpp"
#include "treedata.hpp"
#include "stochasticforest.hpp"
#include "errno.hpp"
#include "options.hpp"
#include "distributions.hpp"

class StochasticForestTest : public CppUnit::TestFixture {
  CPPUNIT_TEST_SUITE( StochasticForestTest );
  CPPUNIT_TEST( test_treeDataPercolation ); 
  CPPUNIT_TEST( test_error );
  CPPUNIT_TEST( test_CART );
  CPPUNIT_TEST_SUITE_END();
  
public:
  
  void setUp();
  void tearDown();

  void test_error();
  void test_treeDataPercolation();

  void test_CART();

private:

  Treedata* trainData_;
  ForestOptions forestOptions_;
  vector<distributions::Random> randoms_;

};

void StochasticForestTest::setUp() {
  
  forestOptions_.setCARTDefaults();
  
  trainData_ = new Treedata("test_103by300_mixed_matrix.afm",'\t',':');

  randoms_.resize(1);
  
}

void StochasticForestTest::tearDown() {

  delete trainData_;
  
}

void StochasticForestTest::test_treeDataPercolation() {


  // First we need some data
  Treedata testData("testdata.tsv",'\t',':');

  // Next load a test predictor
  StochasticForest SF;
  SF.loadForest("test_predictor.sf");

  vector<num_t> prediction,confidence;
  
  SF.predict(&testData,prediction,confidence);

  CPPUNIT_ASSERT( prediction.size() == 4 );
  CPPUNIT_ASSERT( confidence.size() == 4 );
  CPPUNIT_ASSERT( fabs( prediction[0] - 3.9 ) < datadefs::EPS );
  CPPUNIT_ASSERT( fabs( prediction[1] - 4.3 ) < datadefs::EPS );
  CPPUNIT_ASSERT( fabs( prediction[2] - 8.0 ) < datadefs::EPS );
  CPPUNIT_ASSERT( fabs( prediction[3] - 9.0 ) < datadefs::EPS );
  
}


void StochasticForestTest::test_error() {
/*
  options::General_options parameters;
  
  parameters.forestInput = "test_predictor.sf";
  
  // Next load a test predictor
  StochasticForest SF(&parameters);
  
  vector<num_t> x,y;
  CPPUNIT_ASSERT( datadefs::isNAN( SF.error(x,y) ) );
  
  x.push_back(datadefs::NUM_NAN);
  y = x;
  CPPUNIT_ASSERT( datadefs::isNAN( SF.error(x,y) ) );
  
  x.push_back(1.0);
  y.push_back(datadefs::NUM_NAN);
  CPPUNIT_ASSERT( datadefs::isNAN( SF.error(x,y) ) );
  
  x.push_back(0.0);
  y.push_back(2.0);
  CPPUNIT_ASSERT( fabs( SF.error(x,y) - 4.0 ) < datadefs::EPS );
  
  
*/
}


void StochasticForestTest::test_CART() {
 
  Treedata treeData("test_103by300_mixed_nan_matrix.afm",'\t',':');

  forestOptions_.nMaxLeaves = 2;

  StochasticForest CART;

  vector<num_t> featureWeights(treeData.nFeatures(),1);
  featureWeights[0] = 0;

  CART.learnRF(&treeData,0,&forestOptions_,featureWeights,randoms_);

  CPPUNIT_ASSERT( CART.rootNodes_[0]->splitter_.name == "N:input" );

  CPPUNIT_ASSERT( CART.nNodes() == 3 );

  

}


// Registers the fixture into the test 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION( StochasticForestTest ); 

#endif // STOCHASTICFORESTTEST_HPP
