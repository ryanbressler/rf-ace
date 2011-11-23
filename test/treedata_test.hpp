#ifndef TREEDATATEST_HPP
#define TREEDATATEST_HPP

#include <cppunit/extensions/HelperMacros.h>
#include "treedata.hpp"
#include "datadefs.hpp"
#include "errno.hpp"

class treeDataTest : public CppUnit::TestFixture {
  CPPUNIT_TEST_SUITE( treeDataTest );
  CPPUNIT_TEST( test_bracketOperator );
  CPPUNIT_TEST_SUITE_END();
  
public:
  void setUp();
  void tearDown();
  void test_bracketOperator();
};

void treeDataTest::setUp() {}
void treeDataTest::tearDown() {}
 
void treeDataTest::test_bracketOperator() {
  
  string fileName = "test_6by10_featurerows_matrix.tsv";

  Treedata::Treedata treeData(fileName);

  CPPUNIT_ASSERT( treeData.features_.size() == 2*6 );

  for ( size_t i = 0; i < treeData.features_.size(); ++i ) {
    vector<datadefs::num_t> data = treeData.getFeatureData(i);
    for ( size_t j = 0; j < treeData.features_.size(); ++j ) {
      if ( !datadefs::isNAN(data[j]) ) {
	CPPUNIT_ASSERT( fabs(treeData[i][j] - data[j]) < datadefs::EPS );
      }
    }
  }

}


// Registers the fixture into the test 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION( treeDataTest ); 

#endif // TREEDATATEST_HPP
