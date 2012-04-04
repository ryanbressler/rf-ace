#ifndef STOCHASTICFORESTTEST_HPP
#define STOCHASTICFORESTTEST_HPP

#include <cppunit/extensions/HelperMacros.h>
//#include "datadefs.hpp"
//#include "partitionsequence.hpp"
#include "stochasticforest.hpp"
#include "errno.hpp"

class StochasticForestTest : public CppUnit::TestFixture {
  CPPUNIT_TEST_SUITE( StochasticForestTest );
  CPPUNIT_TEST( test_treeDataPercolation ); 
  //CPPUNIT_TEST( test_prediction );
  CPPUNIT_TEST_SUITE_END();
  
public:
  void setUp();
  void tearDown();

  void test_treeDataPercolation();

};

void StochasticForestTest::setUp() {}
void StochasticForestTest::tearDown() {}

void StochasticForestTest::test_treeDataPercolation() {

  // First we need some data
  Treedata::Treedata treeData("testdata.tsv",'\t',':');

  // Next load a test predictor
  StochasticForest SF(&treeData,"test_predictor.sf");

  CPPUNIT_ASSERT( fabs( SF.rootNodes_[0]->getTrainPrediction(0) - 3.9 ) < datadefs::EPS );
  CPPUNIT_ASSERT( fabs( SF.rootNodes_[0]->getTrainPrediction(1) - 4.3 ) < datadefs::EPS );
  CPPUNIT_ASSERT( fabs( SF.rootNodes_[0]->getTrainPrediction(2) - 5.0 ) < datadefs::EPS );
  CPPUNIT_ASSERT( fabs( SF.rootNodes_[0]->getTrainPrediction(3) - 6.6 ) < datadefs::EPS );

  vector<num_t> prediction,confidence;

  SF.predict(prediction,confidence);

  CPPUNIT_ASSERT( prediction.size() == 4 );
  CPPUNIT_ASSERT( confidence.size() == 4 );
  CPPUNIT_ASSERT( fabs( prediction[0] - 3.9 ) < datadefs::EPS );
  CPPUNIT_ASSERT( fabs( prediction[1] - 4.3 ) < datadefs::EPS );
  CPPUNIT_ASSERT( fabs( prediction[2] - 5.0 ) < datadefs::EPS );
  CPPUNIT_ASSERT( fabs( prediction[3] - 6.6 ) < datadefs::EPS );

  SF.rootNodes_[0]->oobIcs_ = utils::range(4);

  prediction = SF.getOobPredictions();

  CPPUNIT_ASSERT( prediction.size() == 4 );
  CPPUNIT_ASSERT( fabs( prediction[0] - 3.9 ) < datadefs::EPS );
  CPPUNIT_ASSERT( fabs( prediction[1] - 4.3 ) < datadefs::EPS );
  CPPUNIT_ASSERT( fabs( prediction[2] - 5.0 ) < datadefs::EPS );
  CPPUNIT_ASSERT( fabs( prediction[3] - 6.6 ) < datadefs::EPS );

}

// Registers the fixture into the test 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION( StochasticForestTest ); 

#endif // STOCHASTICFORESTTEST_HPP
