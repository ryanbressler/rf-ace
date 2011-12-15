#ifndef TREEDATATEST_HPP
#define TREEDATATEST_HPP

#include <cppunit/extensions/HelperMacros.h>
#include "treedata.hpp"
#include "datadefs.hpp"
#include "errno.hpp"

class treeDataTest : public CppUnit::TestFixture {
  CPPUNIT_TEST_SUITE( treeDataTest );
  CPPUNIT_TEST( test_bracketOperator );
  CPPUNIT_TEST( test_name2idxMap );
  CPPUNIT_TEST_SUITE_END();
  
public:
  void setUp();
  void tearDown();
  void test_bracketOperator();
  void test_name2idxMap();
};

void treeDataTest::setUp() {}
void treeDataTest::tearDown() {}
 
void treeDataTest::test_bracketOperator() {
  
  string fileName = "test_6by10_featurerows_matrix.tsv";

  Treedata::Treedata treeData(fileName,'\t',':');

  CPPUNIT_ASSERT( treeData.features_.size() == 2*6 );

  for ( size_t i = 0; i < treeData.features_.size(); ++i ) {
    vector<datadefs::num_t> data = treeData.getFeatureData(i);
    for ( size_t j = 0; j < treeData.features_[0].data.size(); ++j ) {
      if ( !datadefs::isNAN(data[j]) ) {
	CPPUNIT_ASSERT( fabs(treeData[i][j] - data[j]) < datadefs::EPS );
	vector<datadefs::num_t> foo = treeData[treeData.features_[i].name];
	//cout << foo[j] << " vs " << data[j] << endl;
	CPPUNIT_ASSERT( fabs(treeData[treeData.features_[i].name][j] - data[j]) < datadefs::EPS );
      }
    }
  }

}

void treeDataTest::test_name2idxMap() {

  string fileName = "test_6by10_featurerows_matrix.tsv";

  Treedata::Treedata treeData(fileName,'\t',':');

  CPPUNIT_ASSERT( treeData.name2idx_.size() == 2*treeData.nFeatures() );
  CPPUNIT_ASSERT( treeData.features_.size() == 2*treeData.nFeatures() );

  CPPUNIT_ASSERT( treeData.getFeatureName(0) == "N:F1" );
  CPPUNIT_ASSERT( treeData.getFeatureName(1) == "N:F2" );
  CPPUNIT_ASSERT( treeData.getFeatureName(2) == "C:F3" );
  CPPUNIT_ASSERT( treeData.getFeatureName(3) == "N:F4" );
  CPPUNIT_ASSERT( treeData.getFeatureName(4) == "C:F5" );
  CPPUNIT_ASSERT( treeData.getFeatureName(5) == "N:F6" );

  CPPUNIT_ASSERT( treeData.getFeatureName(6) == "N:F1_CONTRAST" );
  CPPUNIT_ASSERT( treeData.getFeatureName(7) == "N:F2_CONTRAST" );
  CPPUNIT_ASSERT( treeData.getFeatureName(8) == "C:F3_CONTRAST" );
  CPPUNIT_ASSERT( treeData.getFeatureName(9) == "N:F4_CONTRAST" );
  CPPUNIT_ASSERT( treeData.getFeatureName(10) == "C:F5_CONTRAST" );
  CPPUNIT_ASSERT( treeData.getFeatureName(11) == "N:F6_CONTRAST" );

  CPPUNIT_ASSERT( treeData.name2idx_["N:F1"] == 0 );
  CPPUNIT_ASSERT( treeData.name2idx_["N:F2"] == 1 );
  CPPUNIT_ASSERT( treeData.name2idx_["C:F3"] == 2 );
  CPPUNIT_ASSERT( treeData.name2idx_["N:F4"] == 3 );
  CPPUNIT_ASSERT( treeData.name2idx_["C:F5"] == 4 );
  CPPUNIT_ASSERT( treeData.name2idx_["N:F6"] == 5 );

  CPPUNIT_ASSERT( treeData.name2idx_["N:F1_CONTRAST"] == 6 );
  CPPUNIT_ASSERT( treeData.name2idx_["N:F2_CONTRAST"] == 7 );
  CPPUNIT_ASSERT( treeData.name2idx_["C:F3_CONTRAST"] == 8 );
  CPPUNIT_ASSERT( treeData.name2idx_["N:F4_CONTRAST"] == 9 );
  CPPUNIT_ASSERT( treeData.name2idx_["C:F5_CONTRAST"] == 10 );
  CPPUNIT_ASSERT( treeData.name2idx_["N:F6_CONTRAST"] == 11 );
  
  

}

// Registers the fixture into the test 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION( treeDataTest ); 

#endif // TREEDATATEST_HPP
