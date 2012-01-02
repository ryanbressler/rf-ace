#ifndef TREEDATATEST_HPP
#define TREEDATATEST_HPP

#include <cppunit/extensions/HelperMacros.h>
#include "treedata.hpp"
#include "datadefs.hpp"
#include "errno.hpp"

class treeDataTest : public CppUnit::TestFixture {
  CPPUNIT_TEST_SUITE( treeDataTest );
  //CPPUNIT_TEST( test_bracketOperator );
  CPPUNIT_TEST( test_name2idxMap );
  CPPUNIT_TEST( test_getFeatureData );
  CPPUNIT_TEST( test_updateSortOrder );
  CPPUNIT_TEST( test_getFilteredAndSortedDataPair2 );
  CPPUNIT_TEST_SUITE_END();
  
public:
  void setUp();
  void tearDown();
  void test_bracketOperator();
  void test_name2idxMap();
  void test_getFeatureData();
  void test_updateSortOrder();
  void test_getFilteredAndSortedDataPair2();
};

void treeDataTest::setUp() {}
void treeDataTest::tearDown() {}
 
/*
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
*/

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

void treeDataTest::test_getFeatureData() {

  Treedata treedata("test_2by8_featurerows_matrix.tsv",'\t',':');

  vector<size_t> sampleIcs(8);
  datadefs::range(sampleIcs);

  vector<num_t> v1,v2;

  treedata.getFilteredDataPair(0,1,sampleIcs,v1,v2);

  

}

void treeDataTest::test_updateSortOrder() {

  Treedata treedata("test_6by10_featurerows_matrix.tsv",'\t',':');
  
  /*
    N:F1    nA      8.5     3.4     7.2     5       6       7       11      9       NA
    N:F2    2       3       4       5       6       NA      NA      9       nan     10
    C:F3    NA      nA      naN     NaN     1       1       1       2       2       2
    N:F4    10      9.9     8       7       6       5       4       3       2.4     1
    C:F5    3       3       3       4       4       5       3       2       2       2
    N:F6    9       8       7       9       8       7       3       2       1.0     99.23
  */

  CPPUNIT_ASSERT( treedata.features_[0].sortOrder[0] == 0 );
  CPPUNIT_ASSERT( treedata.features_[0].sortOrder[1] == 5 );
  CPPUNIT_ASSERT( treedata.features_[0].sortOrder[2] == 0 );
  CPPUNIT_ASSERT( treedata.features_[0].sortOrder[3] == 4 );
  CPPUNIT_ASSERT( treedata.features_[0].sortOrder[4] == 1 );
  CPPUNIT_ASSERT( treedata.features_[0].sortOrder[5] == 2 );
  CPPUNIT_ASSERT( treedata.features_[0].sortOrder[6] == 3 );
  CPPUNIT_ASSERT( treedata.features_[0].sortOrder[7] == 7 );
  CPPUNIT_ASSERT( treedata.features_[0].sortOrder[8] == 6 );
  CPPUNIT_ASSERT( treedata.features_[0].sortOrder[9] == 0 );

  CPPUNIT_ASSERT( treedata.features_[1].sortOrder[0] == 0 );
  CPPUNIT_ASSERT( treedata.features_[1].sortOrder[1] == 1 );
  CPPUNIT_ASSERT( treedata.features_[1].sortOrder[2] == 2 );
  CPPUNIT_ASSERT( treedata.features_[1].sortOrder[3] == 3 );
  CPPUNIT_ASSERT( treedata.features_[1].sortOrder[4] == 4 );
  CPPUNIT_ASSERT( treedata.features_[1].sortOrder[5] == 0 );
  CPPUNIT_ASSERT( treedata.features_[1].sortOrder[6] == 0 );
  CPPUNIT_ASSERT( treedata.features_[1].sortOrder[7] == 5 );
  CPPUNIT_ASSERT( treedata.features_[1].sortOrder[8] == 0 );
  CPPUNIT_ASSERT( treedata.features_[1].sortOrder[9] == 6 );

  CPPUNIT_ASSERT( treedata.features_[2].sortOrder.size() == 0 );

  CPPUNIT_ASSERT( treedata.features_[3].sortOrder[0] == 9 );
  CPPUNIT_ASSERT( treedata.features_[3].sortOrder[1] == 8 );
  CPPUNIT_ASSERT( treedata.features_[3].sortOrder[2] == 7 );
  CPPUNIT_ASSERT( treedata.features_[3].sortOrder[3] == 6 );
  CPPUNIT_ASSERT( treedata.features_[3].sortOrder[4] == 5 );
  CPPUNIT_ASSERT( treedata.features_[3].sortOrder[5] == 4 );
  CPPUNIT_ASSERT( treedata.features_[3].sortOrder[6] == 3 );
  CPPUNIT_ASSERT( treedata.features_[3].sortOrder[7] == 2 );
  CPPUNIT_ASSERT( treedata.features_[3].sortOrder[8] == 1 );
  CPPUNIT_ASSERT( treedata.features_[3].sortOrder[9] == 0 );

  CPPUNIT_ASSERT( treedata.features_[4].sortOrder.size() == 0 );

  CPPUNIT_ASSERT( treedata.features_[5].sortOrder[6] == 2 );
  CPPUNIT_ASSERT( treedata.features_[5].sortOrder[7] == 1 );
  CPPUNIT_ASSERT( treedata.features_[5].sortOrder[8] == 0 );
  CPPUNIT_ASSERT( treedata.features_[5].sortOrder[9] == 9 );

}

void treeDataTest::test_getFilteredAndSortedDataPair2() {

  string fileName = "test_6by10_featurerows_matrix.tsv";

  Treedata::Treedata treeData(fileName,'\t',':');

  vector<size_t> sampleIcs(10);
  vector<num_t> tv,fv;

  datadefs::range(sampleIcs);

  /*
    N:F1    nA      8.5     3.4     7.2     5       6       7       11      9       NA
    N:F2    2       3       4       5       6       NA      NA      9       nan     10
    C:F3    NA      nA      naN     NaN     1       1       1       2       2       2
    N:F4    10      9.9     8       7       6       5       4       3       2.4     1
    C:F5    3       3       3       4       4       5       3       2       2       2
    N:F6    9       8       7       9       8       7       3       2       1.0     99.23
  */

  treeData.getFilteredAndSortedDataPair2(0,1,sampleIcs,tv,fv);

  CPPUNIT_ASSERT( sampleIcs.size() == 5 );
  CPPUNIT_ASSERT( sampleIcs[0] == 1 );
  CPPUNIT_ASSERT( sampleIcs[1] == 2 );
  CPPUNIT_ASSERT( sampleIcs[2] == 3 );
  CPPUNIT_ASSERT( sampleIcs[3] == 4 );
  CPPUNIT_ASSERT( sampleIcs[4] == 7 );

  CPPUNIT_ASSERT( tv.size() == 5 );
  CPPUNIT_ASSERT( fabs( tv[0] - 8.5) < datadefs::EPS );
  CPPUNIT_ASSERT( fabs( tv[1] - 3.4) < datadefs::EPS );
  CPPUNIT_ASSERT( fabs( tv[2] - 7.2) < datadefs::EPS );
  CPPUNIT_ASSERT( fabs( tv[3] - 5) < datadefs::EPS );
  CPPUNIT_ASSERT( fabs( tv[4] - 11) < datadefs::EPS );

  CPPUNIT_ASSERT( fv.size() == 5 );
  CPPUNIT_ASSERT( fabs( fv[0] - 3) < datadefs::EPS );
  CPPUNIT_ASSERT( fabs( fv[1] - 4) < datadefs::EPS );
  CPPUNIT_ASSERT( fabs( fv[2] - 5) < datadefs::EPS );
  CPPUNIT_ASSERT( fabs( fv[3] - 6) < datadefs::EPS );
  CPPUNIT_ASSERT( fabs( fv[4] - 9) < datadefs::EPS );

  sampleIcs.resize(10);
  sampleIcs[0] = 3; // 5
  sampleIcs[1] = 4; // 3
  sampleIcs[2] = 3; // 5
  sampleIcs[3] = 5; // 4
  sampleIcs[4] = 0; // NA
  sampleIcs[5] = 8; // 6
  sampleIcs[6] = 2; // 2
  sampleIcs[7] = 2; // 2
  sampleIcs[8] = 3; // 5 
  sampleIcs[9] = 0; // NA

  // multiplicity: 3.4(2) 5(1) 6(1) 7(0) 7.2(3) 8.5(0) 9(1) 11(0) 0 0 

  /*
    N:F1    nA      8.5     3.4     7.2     5       6       7       11      9       NA
    N:F2    2       3       4       5       6       NA      NA      9       nan     10
    C:F3    NA      nA      naN     NaN     1       1       1       2       2       2
    N:F4    10      9.9     8       7       6       5       4       3       2.4     1
    C:F5    3       3       3       4       4       5       3       2       2       2
    N:F6    9       8       7       9       8       7       3       2       1.0     99.23
  */

  treeData.getFilteredAndSortedDataPair2(5,0,sampleIcs,tv,fv);

  //datadefs::print(sampleIcs);

  CPPUNIT_ASSERT( sampleIcs.size() ==  8 );
  CPPUNIT_ASSERT( sampleIcs[0] == 2 );
  CPPUNIT_ASSERT( sampleIcs[1] == 2 );
  CPPUNIT_ASSERT( sampleIcs[2] == 4 );
  CPPUNIT_ASSERT( sampleIcs[3] == 5 );
  CPPUNIT_ASSERT( sampleIcs[4] == 3 );
  CPPUNIT_ASSERT( sampleIcs[5] == 3 );
  CPPUNIT_ASSERT( sampleIcs[6] == 3 );
  CPPUNIT_ASSERT( sampleIcs[7] == 8 );

  CPPUNIT_ASSERT( tv.size() == 8 );
  CPPUNIT_ASSERT( fabs( tv[0] - 7 ) < datadefs::EPS );
  CPPUNIT_ASSERT( fabs( tv[1] - 7 ) < datadefs::EPS );
  CPPUNIT_ASSERT( fabs( tv[2] - 8 ) < datadefs::EPS );
  CPPUNIT_ASSERT( fabs( tv[3] - 7 ) < datadefs::EPS );
  CPPUNIT_ASSERT( fabs( tv[4] - 9 ) < datadefs::EPS );
  CPPUNIT_ASSERT( fabs( tv[5] - 9 ) < datadefs::EPS );
  CPPUNIT_ASSERT( fabs( tv[6] - 9 ) < datadefs::EPS );
  CPPUNIT_ASSERT( fabs( tv[7] - 1.0 ) < datadefs::EPS );

  CPPUNIT_ASSERT( fv.size() == 8 );
  CPPUNIT_ASSERT( fabs( fv[0] - 3.4 ) < datadefs::EPS );
  CPPUNIT_ASSERT( fabs( fv[1] - 3.4 ) < datadefs::EPS );
  CPPUNIT_ASSERT( fabs( fv[2] - 5 ) < datadefs::EPS );
  CPPUNIT_ASSERT( fabs( fv[3] - 6 ) < datadefs::EPS );
  CPPUNIT_ASSERT( fabs( fv[4] - 7.2 ) < datadefs::EPS );
  CPPUNIT_ASSERT( fabs( fv[5] - 7.2 ) < datadefs::EPS );
  CPPUNIT_ASSERT( fabs( fv[6] - 7.2 ) < datadefs::EPS );
  CPPUNIT_ASSERT( fabs( fv[7] - 9 ) < datadefs::EPS );


}


// Registers the fixture into the test 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION( treeDataTest ); 

#endif // TREEDATATEST_HPP
