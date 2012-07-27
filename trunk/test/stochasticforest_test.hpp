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

  //Treedata trainData_;
  //options::General_options* genOp_;
  //StochasticForest* CART_;

};

void StochasticForestTest::setUp() {

  /*
    trainData_ = new Treedata("test_103by300_mixed_matrix.afm",'\t',':');
    
    genOp_ = new options::General_options;
    
    genOp_->setCARTDefaults();
    
    genOp_->nTrees = 1;
    genOp_->nMaxLeaves = 1000;
    genOp_->nodeSize = 5;
    
    CART_ = new StochasticForest(trainData_,genOp_);
  */

}

void StochasticForestTest::tearDown() {

  //delete CART_;
  //delete trainData_;
  //delete genOp_;
  
}

void StochasticForestTest::test_treeDataPercolation() {


  options::General_options parameters;
  //parameters.modelType = options::RF;
  parameters.loadUserParams();
  parameters.forestInput = "test_predictor.sf";

  // First we need some data
  Treedata::Treedata testData("testdata.tsv",'\t',':',parameters.randIntGens[0]);
  
  // Next load a test predictor
  StochasticForest SF(&parameters);

  vector<num_t> prediction,confidence;
  
  SF.predictWithTestData(&testData,prediction,confidence);
  
  CPPUNIT_ASSERT( prediction.size() == 4 );
  CPPUNIT_ASSERT( confidence.size() == 4 );
  CPPUNIT_ASSERT( fabs( prediction[0] - 3.9 ) < datadefs::EPS );
  CPPUNIT_ASSERT( fabs( prediction[1] - 4.3 ) < datadefs::EPS );
  CPPUNIT_ASSERT( fabs( prediction[2] - 3.9 ) < datadefs::EPS );
  CPPUNIT_ASSERT( fabs( prediction[3] - 6.5 ) < datadefs::EPS );
  

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
 
  
  options::General_options CARTParameters;
  CARTParameters.modelType = options::CART;
  CARTParameters.loadUserParams();
  CARTParameters.targetStr = "N:output";

  Treedata treeData("test_103by300_mixed_matrix.afm",'\t',':',CARTParameters.randIntGens[0]);

  CARTParameters.nMaxLeaves = treeData.nSamples();

  //CARTParameters.printParameters();

  StochasticForest CART(&treeData,&CARTParameters);

}


// Registers the fixture into the test 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION( StochasticForestTest ); 

#endif // STOCHASTICFORESTTEST_HPP
