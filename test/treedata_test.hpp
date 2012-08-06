#ifndef TREEDATATEST_HPP
#define TREEDATATEST_HPP

#include <cppunit/extensions/HelperMacros.h>
#include "treedata.hpp"
#include "datadefs.hpp"
#include "node.hpp"
#include "errno.hpp"
#include "distributions.hpp"

class treeDataTest : public CppUnit::TestFixture {
  CPPUNIT_TEST_SUITE( treeDataTest );
  CPPUNIT_TEST( test_name2idxMap );
  CPPUNIT_TEST( test_getFeatureData );
  CPPUNIT_TEST( test_getFilteredFeatureData );
  CPPUNIT_TEST( test_getFilteredAndSortedFeatureDataPair3 );
  CPPUNIT_TEST( test_parseARFF );
  CPPUNIT_TEST( test_parseARFF_extended ); 
  CPPUNIT_TEST( test_keepFeatures );
  CPPUNIT_TEST( test_removeFeatures );
  CPPUNIT_TEST( test_numericalFeatureSplitsNumericalTarget );
  CPPUNIT_TEST( test_numericalFeatureSplitsCategoricalTarget );
  CPPUNIT_TEST( test_categoricalFeatureSplitsNumericalTarget );
  CPPUNIT_TEST( test_categoricalFeatureSplitsCategoricalTarget );
  CPPUNIT_TEST( test_fullSplitterSweep );
  CPPUNIT_TEST( test_end );
  CPPUNIT_TEST_SUITE_END();
  
public:
  void setUp();
  void tearDown();
  void test_permuteContrasts();
  void test_name2idxMap();
  void test_getFeatureData();
  void test_getFilteredFeatureData();
  void test_getFilteredAndSortedFeatureDataPair3();
  void test_parseARFF();
  void test_parseARFF_extended();
  void test_keepFeatures();
  void test_removeFeatures();
  void test_numericalFeatureSplitsNumericalTarget();
  void test_numericalFeatureSplitsCategoricalTarget();
  void test_categoricalFeatureSplitsNumericalTarget();
  void test_categoricalFeatureSplitsCategoricalTarget();
  void test_fullSplitterSweep();
  void test_end();

private:

  Treedata* treeData_;
  Treedata* treeData_103by300NaN_;

  distributions::RandInt randInt_;



};

void treeDataTest::setUp() {

  randInt_.seed(0);

  treeData_ = new Treedata("test_103by300_mixed_matrix.afm",'\t',':',randInt_);
  treeData_103by300NaN_ = new Treedata("test_103by300_mixed_nan_matrix.afm",'\t',':',randInt_);

}

void treeDataTest::tearDown() {
  
  delete treeData_;
  delete treeData_103by300NaN_;
  
}
 

void treeDataTest::test_permuteContrasts() {
  
  string fileName = "test_6by10_mixed_matrix.tsv";
  
  Treedata::Treedata treeData(fileName,'\t',':',randInt_);
  
  //treeData.permuteContrasts();

  CPPUNIT_ASSERT( treeData.features_.size() == 2*6 );
  
  /*
    N:F1    nA      8.5     3.4     7.2     5       6       7       11      9       NA
    N:F2    2       3       4       5       6       NA      NA      9       nan     10
    C:F3    NA      nA      naN     NaN     1       1       1       2       2       2
    N:F4    10      9.9     8       7       6       5       4       3       2.4     1
    C:F5    3       3       3       4       4       5       3       2       2       2
    N:F6    9       8       7       9       8       7       3       2       1.0     99.23
  */
  
  CPPUNIT_ASSERT( treeData.nRealSamples(9) == 10 );
  CPPUNIT_ASSERT( treeData.nRealSamples(10) == 10 );
  CPPUNIT_ASSERT( treeData.nRealSamples(11) == 10 );

  vector<num_t> featureData = treeData.getFeatureData(6);
  CPPUNIT_ASSERT( treeData.nRealSamples(6) == 8 );
  CPPUNIT_ASSERT( datadefs::isNAN(featureData[0]) );
  CPPUNIT_ASSERT( datadefs::isNAN(featureData[9]) );

  featureData = treeData.getFeatureData(7);
  CPPUNIT_ASSERT( treeData.nRealSamples(7) == 7 );
  CPPUNIT_ASSERT( datadefs::isNAN(featureData[5]) );
  CPPUNIT_ASSERT( datadefs::isNAN(featureData[6]) );
  CPPUNIT_ASSERT( datadefs::isNAN(featureData[8]) );

  featureData = treeData.getFeatureData(8);
  CPPUNIT_ASSERT( treeData.nRealSamples(8) == 6 );
  CPPUNIT_ASSERT( datadefs::isNAN(featureData[0]) );
  CPPUNIT_ASSERT( datadefs::isNAN(featureData[1]) );
  CPPUNIT_ASSERT( datadefs::isNAN(featureData[2]) );
  CPPUNIT_ASSERT( datadefs::isNAN(featureData[3]) );

}

void treeDataTest::test_name2idxMap() {

  string fileName = "test_6by10_mixed_matrix.tsv";

  Treedata::Treedata treeData(fileName,'\t',':',randInt_);

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

  Treedata treedata("test_2by8_numerical_matrix.tsv",'\t',':',randInt_);

  vector<size_t> sampleIcs = utils::range(8);

  vector<num_t> v1,v2;

  treedata.getFilteredFeatureDataPair(0,1,sampleIcs,v1,v2);

  

}


void treeDataTest::test_getFilteredFeatureData() {
  
  string fileName = "test_6by10_mixed_matrix.tsv";
  
  Treedata::Treedata treeData(fileName,'\t',':',randInt_);


  /*
    N:F1    nA      8.5     3.4     7.2     5       6       7       11      9       NA
    N:F2    2       3       4       5       6       NA      NA      9       nan     10
    C:F3    NA      nA      naN     NaN     1       1       1       2       2       2
    N:F4    10      9.9     8       7       6       5       4       3       2.4     1
    C:F5    3       3       3       4       4       5       3       2       2       2
    N:F6    9       8       7       9       8       7       3       2       1.0     99.23
  */
  
  vector<size_t> sampleIcs = utils::range(10);

  vector<num_t> filteredData = treeData.getFilteredFeatureData(0,sampleIcs);

  CPPUNIT_ASSERT( filteredData.size() == 8 );
  CPPUNIT_ASSERT( fabs( filteredData[0] - 8.5 ) < datadefs::EPS );
  CPPUNIT_ASSERT( fabs( filteredData[1] - 3.4 ) < datadefs::EPS );
  CPPUNIT_ASSERT( fabs( filteredData[2] - 7.2 ) < datadefs::EPS );
  CPPUNIT_ASSERT( fabs( filteredData[3] - 5 )   < datadefs::EPS );
  CPPUNIT_ASSERT( fabs( filteredData[4] - 6 )   < datadefs::EPS );
  CPPUNIT_ASSERT( fabs( filteredData[5] - 7 )   < datadefs::EPS );
  CPPUNIT_ASSERT( fabs( filteredData[6] - 11 )  < datadefs::EPS );
  CPPUNIT_ASSERT( fabs( filteredData[7] - 9 )   < datadefs::EPS );

  CPPUNIT_ASSERT( sampleIcs.size() == 8 );
  CPPUNIT_ASSERT( sampleIcs[0] == 1 );
  CPPUNIT_ASSERT( sampleIcs[1] == 2 );
  CPPUNIT_ASSERT( sampleIcs[2] == 3 );
  CPPUNIT_ASSERT( sampleIcs[3] == 4 );
  CPPUNIT_ASSERT( sampleIcs[4] == 5 );
  CPPUNIT_ASSERT( sampleIcs[5] == 6 );
  CPPUNIT_ASSERT( sampleIcs[6] == 7 );
  CPPUNIT_ASSERT( sampleIcs[7] == 8 );

  filteredData = treeData.getFilteredFeatureData(1,sampleIcs);

  CPPUNIT_ASSERT( filteredData.size() == 5 );
  CPPUNIT_ASSERT( fabs( filteredData[0] - 3 ) < datadefs::EPS );
  CPPUNIT_ASSERT( fabs( filteredData[1] - 4 ) < datadefs::EPS );
  CPPUNIT_ASSERT( fabs( filteredData[2] - 5 ) < datadefs::EPS );
  CPPUNIT_ASSERT( fabs( filteredData[3] - 6 ) < datadefs::EPS );
  CPPUNIT_ASSERT( fabs( filteredData[4] - 9 ) < datadefs::EPS );

  CPPUNIT_ASSERT( sampleIcs.size() == 5 );
  CPPUNIT_ASSERT( sampleIcs[0] == 1 );
  CPPUNIT_ASSERT( sampleIcs[1] == 2 );
  CPPUNIT_ASSERT( sampleIcs[2] == 3 );
  CPPUNIT_ASSERT( sampleIcs[3] == 4 );
  CPPUNIT_ASSERT( sampleIcs[4] == 7 );

  sampleIcs[0] = 0;
  sampleIcs[1] = 0;
  sampleIcs[2] = 0;
  sampleIcs[3] = 5;
  sampleIcs[4] = 5;

  filteredData = treeData.getFilteredFeatureData(0,sampleIcs);

  CPPUNIT_ASSERT( filteredData.size() == 2 );
  CPPUNIT_ASSERT( fabs( filteredData[0] - 6 ) < datadefs::EPS );
  CPPUNIT_ASSERT( fabs( filteredData[1] - 6 ) < datadefs::EPS );

  CPPUNIT_ASSERT( sampleIcs.size() == 2 );
  CPPUNIT_ASSERT( sampleIcs[0] == 5 );
  CPPUNIT_ASSERT( sampleIcs[1] == 5 );

  filteredData = treeData.getFilteredFeatureData(1,sampleIcs);

  CPPUNIT_ASSERT( filteredData.size() == 0 );
  CPPUNIT_ASSERT( sampleIcs.size() == 0 );
  
 


}



void treeDataTest::test_getFilteredAndSortedFeatureDataPair3() {

  string fileName = "test_6by10_mixed_matrix.tsv";

  Treedata::Treedata treeData(fileName,'\t',':',randInt_);

  /*
    N:F1    nA      8.5     3.4     7.2     5       6       7       11      9       NA
    N:F2    2       3       4       5       6       NA      NA      9       nan     10
    C:F3    NA      nA      naN     NaN     1       1       1       2       2       2
    N:F4    10      9.9     8       7       6       5       4       3       2.4     1
    C:F5    3       3       3       4       4       5       3       2       2       2
    N:F6    9       8       7       9       8       7       3       2       1.0     99.23
  */


  vector<size_t> sampleIcs;
  sampleIcs.push_back(1);
  sampleIcs.push_back(2);
  sampleIcs.push_back(3);
  sampleIcs.push_back(4);
  sampleIcs.push_back(5);
  sampleIcs.push_back(6);
  sampleIcs.push_back(7);
  sampleIcs.push_back(8);

  vector<num_t> tv,fv;

  //datadefs::range(sampleIcs);

  treeData.getFilteredAndSortedFeatureDataPair3(0,1,sampleIcs,tv,fv);

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

  treeData.getFilteredAndSortedFeatureDataPair3(5,0,sampleIcs,tv,fv);

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

void treeDataTest::test_parseARFF() {
  
  Treedata treeData("test_5by10_numeric_matrix.arff",'\t',':',randInt_);

  /*
    0.8147, 1.0000,  0.0596, 0.9160, 6.0000
    0.9058, 2.0000,  0.6820, 0.0012, 14.0000
    0.1270, 3.0000,  0.0424, 0.4624, 24.0000
    0.9134, 4.0000,  0.0714, 0.4243, 36.0000
    0.6324, 5.0000,  ?,      0.4609, 50.0000
    0.0975, 6.0000,  0.0967, 0.7702, 66.0000
    0.2785, 7.0000,  0.8181, 0.3225, 84.0000
    0.5469, ?,       0.8175, 0.7847, 104.0000
    0.9575, 9.0000,  0.7224, 0.4714, 126.0000
    0.9649, 10.0000, 0.1499, 0.0358, 150.0000
  */

  CPPUNIT_ASSERT( treeData.nFeatures() == 5 );
  CPPUNIT_ASSERT( treeData.nSamples() == 10 );

  CPPUNIT_ASSERT( treeData.nRealSamples(0) == 10 );
  CPPUNIT_ASSERT( treeData.nRealSamples(1) == 9 );
  CPPUNIT_ASSERT( treeData.nRealSamples(2) == 9 );
  CPPUNIT_ASSERT( treeData.nRealSamples(3) == 10 );
  CPPUNIT_ASSERT( treeData.nRealSamples(4) == 10 );

  CPPUNIT_ASSERT( treeData.nRealSamples(0,1) == 9 );
  CPPUNIT_ASSERT( treeData.nRealSamples(1,2) == 8 );

  CPPUNIT_ASSERT( treeData.getFeatureName(0) == "x1" );
  CPPUNIT_ASSERT( treeData.getFeatureName(1) == "x2" );
  CPPUNIT_ASSERT( treeData.getFeatureName(2) == "x3" );
  CPPUNIT_ASSERT( treeData.getFeatureName(3) == "x4" );
  CPPUNIT_ASSERT( treeData.getFeatureName(4) == "y"  );

  CPPUNIT_ASSERT( fabs( treeData.getFeatureData(0,0) - 0.8147 ) < datadefs::EPS );
  CPPUNIT_ASSERT( fabs( treeData.getFeatureData(0,1) - 0.9058 ) < datadefs::EPS );
  CPPUNIT_ASSERT( fabs( treeData.getFeatureData(0,2) - 0.1270 ) < datadefs::EPS );
  CPPUNIT_ASSERT( fabs( treeData.getFeatureData(0,3) - 0.9134 ) < datadefs::EPS );
  CPPUNIT_ASSERT( fabs( treeData.getFeatureData(0,4) - 0.6324 ) < datadefs::EPS );
  CPPUNIT_ASSERT( fabs( treeData.getFeatureData(0,5) - 0.0975 ) < datadefs::EPS );
  CPPUNIT_ASSERT( fabs( treeData.getFeatureData(0,6) - 0.2785 ) < datadefs::EPS );
  CPPUNIT_ASSERT( fabs( treeData.getFeatureData(0,7) - 0.5469 ) < datadefs::EPS );
  CPPUNIT_ASSERT( fabs( treeData.getFeatureData(0,8) - 0.9575 ) < datadefs::EPS );
  CPPUNIT_ASSERT( fabs( treeData.getFeatureData(0,9) - 0.9649 ) < datadefs::EPS );

  CPPUNIT_ASSERT( fabs( treeData.getFeatureData(1,0) - 1.0000 ) < datadefs::EPS );
  CPPUNIT_ASSERT( fabs( treeData.getFeatureData(1,1) - 2.0000 ) < datadefs::EPS );
  CPPUNIT_ASSERT( fabs( treeData.getFeatureData(1,2) - 3.0000 ) < datadefs::EPS );
  CPPUNIT_ASSERT( fabs( treeData.getFeatureData(1,3) - 4.0000 ) < datadefs::EPS );
  CPPUNIT_ASSERT( fabs( treeData.getFeatureData(1,4) - 5.0000 ) < datadefs::EPS );
  CPPUNIT_ASSERT( fabs( treeData.getFeatureData(1,5) - 6.0000 ) < datadefs::EPS );
  CPPUNIT_ASSERT( fabs( treeData.getFeatureData(1,6) - 7.0000 ) < datadefs::EPS );
  CPPUNIT_ASSERT( datadefs::isNAN(treeData.getFeatureData(1,7)) );
  CPPUNIT_ASSERT( fabs( treeData.getFeatureData(1,8) - 9.0000 ) < datadefs::EPS );
  CPPUNIT_ASSERT( fabs( treeData.getFeatureData(1,9) - 10.0000 ) < datadefs::EPS );

  /*
    0.8147, 1.0000,  0.0596, 0.9160, 6.0000
    0.9058, 2.0000,  0.6820, 0.0012, 14.0000
    0.1270, 3.0000,  0.0424, 0.4624, 24.0000
    0.9134, 4.0000,  0.0714, 0.4243, 36.0000
    0.6324, 5.0000,  ?,      0.4609, 50.0000
    0.0975, 6.0000,  0.0967, 0.7702, 66.0000
    0.2785, 7.0000,  0.8181, 0.3225, 84.0000
    0.5469, ?,       0.8175, 0.7847, 104.0000
    0.9575, 9.0000,  0.7224, 0.4714, 126.0000
    0.9649, 10.0000, 0.1499, 0.0358, 150.0000
  */

  CPPUNIT_ASSERT( fabs( treeData.getFeatureData(2,0) - 0.0596 ) < datadefs::EPS );
  CPPUNIT_ASSERT( fabs( treeData.getFeatureData(2,1) - 0.6820 ) < datadefs::EPS );
  CPPUNIT_ASSERT( fabs( treeData.getFeatureData(2,2) - 0.0424 ) < datadefs::EPS );
  CPPUNIT_ASSERT( fabs( treeData.getFeatureData(2,3) - 0.0714 ) < datadefs::EPS );
  CPPUNIT_ASSERT( datadefs::isNAN(treeData.getFeatureData(2,4)) );
  CPPUNIT_ASSERT( fabs( treeData.getFeatureData(2,5) - 0.0967 ) < datadefs::EPS );
  CPPUNIT_ASSERT( fabs( treeData.getFeatureData(2,6) - 0.8181 ) < datadefs::EPS );
  CPPUNIT_ASSERT( fabs( treeData.getFeatureData(2,7) - 0.8175 ) < datadefs::EPS );
  CPPUNIT_ASSERT( fabs( treeData.getFeatureData(2,8) - 0.7224 ) < datadefs::EPS );
  CPPUNIT_ASSERT( fabs( treeData.getFeatureData(2,9) - 0.1499 ) < datadefs::EPS );

  CPPUNIT_ASSERT( fabs( treeData.getFeatureData(3,0) - 0.9160 ) < datadefs::EPS );
  CPPUNIT_ASSERT( fabs( treeData.getFeatureData(3,1) - 0.0012 ) < datadefs::EPS );
  CPPUNIT_ASSERT( fabs( treeData.getFeatureData(3,2) - 0.4624 ) < datadefs::EPS );
  CPPUNIT_ASSERT( fabs( treeData.getFeatureData(3,3) - 0.4243 ) < datadefs::EPS );
  CPPUNIT_ASSERT( fabs( treeData.getFeatureData(3,4) - 0.4609 ) < datadefs::EPS );
  CPPUNIT_ASSERT( fabs( treeData.getFeatureData(3,5) - 0.7702 ) < datadefs::EPS );
  CPPUNIT_ASSERT( fabs( treeData.getFeatureData(3,6) - 0.3225 ) < datadefs::EPS );
  CPPUNIT_ASSERT( fabs( treeData.getFeatureData(3,7) - 0.7847 ) < datadefs::EPS );
  CPPUNIT_ASSERT( fabs( treeData.getFeatureData(3,8) - 0.4714 ) < datadefs::EPS );
  CPPUNIT_ASSERT( fabs( treeData.getFeatureData(3,9) - 0.0358 ) < datadefs::EPS );

  CPPUNIT_ASSERT( fabs( treeData.getFeatureData(4,0) - 6.0000 ) < datadefs::EPS );
  CPPUNIT_ASSERT( fabs( treeData.getFeatureData(4,1) - 14.0000 ) < datadefs::EPS );
  CPPUNIT_ASSERT( fabs( treeData.getFeatureData(4,2) - 24.0000 ) < datadefs::EPS );
  CPPUNIT_ASSERT( fabs( treeData.getFeatureData(4,3) - 36.0000 ) < datadefs::EPS );
  CPPUNIT_ASSERT( fabs( treeData.getFeatureData(4,4) - 50.0000 ) < datadefs::EPS );
  CPPUNIT_ASSERT( fabs( treeData.getFeatureData(4,5) - 66.0000 ) < datadefs::EPS );
  CPPUNIT_ASSERT( fabs( treeData.getFeatureData(4,6) - 84.0000 ) < datadefs::EPS );
  CPPUNIT_ASSERT( fabs( treeData.getFeatureData(4,7) - 104.0000 ) < datadefs::EPS );
  CPPUNIT_ASSERT( fabs( treeData.getFeatureData(4,8) - 126.0000 ) < datadefs::EPS );
  CPPUNIT_ASSERT( fabs( treeData.getFeatureData(4,9) - 150.0000 ) < datadefs::EPS );

}

void treeDataTest::test_parseARFF_extended() {

  Treedata treeData("test_12by21_categorical_matrix.arff",'\t',':',randInt_);

}

void treeDataTest::test_keepFeatures() {
  
  Treedata treedata("test_6by10_mixed_matrix.tsv",'\t',':',randInt_);

  // Feature names in the matrix are N:F1 N:F2 C:F3 N:F4 C:F5 N:F6

  set<string> keepFeatureNames;
  keepFeatureNames.insert("N:F2");
  keepFeatureNames.insert("C:F5");

  treedata.whiteList(keepFeatureNames);

  size_t featureIdx;
  featureIdx = treedata.getFeatureIdx("N:F2");
  featureIdx = treedata.getFeatureIdx("C:F5");

  CPPUNIT_ASSERT( treedata.nFeatures() == keepFeatureNames.size() );
  CPPUNIT_ASSERT( treedata.nSamples() == 10 );

  // Internally, treedata also stores the contrasts, doubling the number of features
  CPPUNIT_ASSERT( treedata.name2idx_.size() == 2*treedata.nFeatures() );

  for ( size_t i = 0; i < treedata.nFeatures(); ++i ) {
    string featureName = treedata.getFeatureName(i);
    CPPUNIT_ASSERT( featureName.append("_CONTRAST") == treedata.getFeatureName( treedata.nFeatures() + i ));
  }
 

}

void treeDataTest::test_removeFeatures() {

  Treedata treedata("test_6by10_mixed_matrix.tsv",'\t',':',randInt_);

  // Feature names in the matrix are N:F1 N:F2 C:F3 N:F4 C:F5 N:F6

  set<string> removeFeatureNames;
  removeFeatureNames.insert("N:F2");
  removeFeatureNames.insert("C:F5");
  
  size_t nFeaturesOld = treedata.nFeatures();

  treedata.blackList(removeFeatureNames);

  size_t featureIdx;
  featureIdx = treedata.getFeatureIdx("N:F1");
  featureIdx = treedata.getFeatureIdx("C:F3");
  featureIdx = treedata.getFeatureIdx("N:F4");
  featureIdx = treedata.getFeatureIdx("N:F6");

  CPPUNIT_ASSERT( treedata.nFeatures() == nFeaturesOld - removeFeatureNames.size() );
  CPPUNIT_ASSERT( treedata.nSamples() == 10 );

  // Internally, treedata also stores the contrasts, doubling the number of features
  CPPUNIT_ASSERT( treedata.name2idx_.size() == 2*treedata.nFeatures() );

  for ( size_t i = 0; i < treedata.nFeatures(); ++i ) {
    string featureName = treedata.getFeatureName(i);
    CPPUNIT_ASSERT( featureName.append("_CONTRAST") == treedata.getFeatureName( treedata.nFeatures() + i ));
  }


}

void treeDataTest::test_numericalFeatureSplitsNumericalTarget() {

  vector<size_t> sampleIcs_left(0);
  vector<size_t> sampleIcs_right = utils::range(300);

  datadefs::num_t splitValue;
  datadefs::num_t deltaImpurity;

  size_t minSamples = 1;

  size_t targetIdx = 0; // numerical
  size_t featureIdx = 2; // numerical

  deltaImpurity = treeData_->numericalFeatureSplit(targetIdx,
						   featureIdx,
						   minSamples,
						   sampleIcs_left,
						   sampleIcs_right,
						   splitValue);
  
  {
    set<size_t> leftIcs(sampleIcs_left.begin(),sampleIcs_left.end());
    set<size_t> rightIcs(sampleIcs_right.begin(),sampleIcs_right.end());
    
    CPPUNIT_ASSERT( fabs( deltaImpurity - 1.289680394982406 ) < 1e-10 );
    CPPUNIT_ASSERT( fabs( splitValue - 4.387 ) < 1e-10 );
    
    CPPUNIT_ASSERT( sampleIcs_left.size() == 127 );
    CPPUNIT_ASSERT( sampleIcs_right.size() == 173 );

    CPPUNIT_ASSERT( leftIcs.find(198) != leftIcs.end() );
    CPPUNIT_ASSERT( leftIcs.find(8)   != leftIcs.end() );
    CPPUNIT_ASSERT( leftIcs.find(102) != leftIcs.end() );
    CPPUNIT_ASSERT( leftIcs.find(4)   != leftIcs.end() );
    CPPUNIT_ASSERT( leftIcs.find(299) != leftIcs.end() );
    
    CPPUNIT_ASSERT( rightIcs.find(26) != rightIcs.end() );
    CPPUNIT_ASSERT( rightIcs.find(2)  != rightIcs.end() );
    CPPUNIT_ASSERT( rightIcs.find(10) != rightIcs.end() );
    CPPUNIT_ASSERT( rightIcs.find(81) != rightIcs.end() );
    CPPUNIT_ASSERT( rightIcs.find(33) != rightIcs.end() );
  }
  
  minSamples = 50;

  

}

void treeDataTest::test_numericalFeatureSplitsCategoricalTarget() {

  size_t targetIdx = 1; // categorical
  size_t featureIdx = 2; // numerical

  vector<size_t> sampleIcs_left;
  vector<size_t> sampleIcs_right = utils::range(300);

  datadefs::num_t splitValue;
  datadefs::num_t deltaImpurity;

  size_t minSamples = 1;

  deltaImpurity = treeData_->numericalFeatureSplit(targetIdx,
                                                   featureIdx,
                                                   minSamples,
                                                   sampleIcs_left,
                                                   sampleIcs_right,
                                                   splitValue);

  {
    set<size_t> leftIcs(sampleIcs_left.begin(),sampleIcs_left.end());
    set<size_t> rightIcs(sampleIcs_right.begin(),sampleIcs_right.end());

    CPPUNIT_ASSERT( fabs( deltaImpurity - 0.012389077212806 ) < 1e-10 );
    CPPUNIT_ASSERT( fabs( splitValue - 9.827 ) < 1e-10 );

    CPPUNIT_ASSERT( sampleIcs_left.size() == 295 );
    CPPUNIT_ASSERT( sampleIcs_right.size() == 5 );

    CPPUNIT_ASSERT( leftIcs.find(261) != leftIcs.end() );
    CPPUNIT_ASSERT( leftIcs.find(185) != leftIcs.end() );
    CPPUNIT_ASSERT( leftIcs.find(3)   != leftIcs.end() );
    CPPUNIT_ASSERT( leftIcs.find(7)   != leftIcs.end() );
    CPPUNIT_ASSERT( leftIcs.find(256) != leftIcs.end() );

    CPPUNIT_ASSERT( rightIcs.find(69)  != rightIcs.end() );
    CPPUNIT_ASSERT( rightIcs.find(55)  != rightIcs.end() );
    CPPUNIT_ASSERT( rightIcs.find(100) != rightIcs.end() );
    CPPUNIT_ASSERT( rightIcs.find(127) != rightIcs.end() );
    CPPUNIT_ASSERT( rightIcs.find(91)  != rightIcs.end() );
  }


}

void treeDataTest::test_categoricalFeatureSplitsNumericalTarget() {

  vector<size_t> sampleIcs_left(0);
  vector<size_t> sampleIcs_right = utils::range(300);

  set<num_t> splitValues_left,splitValues_right;
  datadefs::num_t deltaImpurity;

  size_t featureIdx = 1;
  size_t targetIdx = 0;
  size_t minSamples = 1;

  deltaImpurity = treeData_->categoricalFeatureSplit(targetIdx,
						     featureIdx,
						     minSamples,
						     sampleIcs_left,
						     sampleIcs_right,
						     splitValues_left,
						     splitValues_right);

  set<string> rawSplitValues_left,rawSplitValues_right;

  for(set<num_t>::const_iterator it(splitValues_left.begin()); it != splitValues_left.end(); ++it ) {
    rawSplitValues_left.insert(treeData_->getRawFeatureData(featureIdx,*it));
  }

  for(set<num_t>::const_iterator it(splitValues_right.begin()); it != splitValues_right.end(); ++it ) {
    rawSplitValues_right.insert(treeData_->getRawFeatureData(featureIdx,*it));
  }

  datadefs::num_t leftFraction = 1.0*sampleIcs_left.size() / ( sampleIcs_left.size() + sampleIcs_right.size() );

  Node node;
  node.setSplitter(0,"foo",leftFraction,rawSplitValues_left,rawSplitValues_right);

  CPPUNIT_ASSERT( node.percolateData("1") == node.rightChild() );
  CPPUNIT_ASSERT( node.percolateData("2") == node.rightChild() );
  CPPUNIT_ASSERT( node.percolateData("3") == node.leftChild() );

  CPPUNIT_ASSERT( sampleIcs_left.size() == 91 );
  CPPUNIT_ASSERT( sampleIcs_right.size() == 209 );

  CPPUNIT_ASSERT( fabs( deltaImpurity - 1.102087375288799 ) < 1e-10 );

}

void treeDataTest::test_categoricalFeatureSplitsCategoricalTarget() {

  vector<size_t> sampleIcs_left(0);
  vector<size_t> sampleIcs_right = utils::range(300);

  set<num_t> splitValues_left,splitValues_right;
  datadefs::num_t deltaImpurity;

  size_t featureIdx = 8;
  size_t targetIdx = 1;
  size_t minSamples = 1;

  deltaImpurity = treeData_->categoricalFeatureSplit(targetIdx,
                                                     featureIdx,
                                                     minSamples,
                                                     sampleIcs_left,
                                                     sampleIcs_right,
                                                     splitValues_left,
                                                     splitValues_right);

  CPPUNIT_ASSERT( sampleIcs_left.size() == 89 );
  CPPUNIT_ASSERT( sampleIcs_right.size() == 211 );

  set<string> rawSplitValues_left,rawSplitValues_right;
  
  for(set<num_t>::const_iterator it(splitValues_left.begin()); it != splitValues_left.end(); ++it ) {
    rawSplitValues_left.insert(treeData_->getRawFeatureData(featureIdx,*it));
  }

  for(set<num_t>::const_iterator it(splitValues_right.begin()); it != splitValues_right.end(); ++it ) {
    rawSplitValues_right.insert(treeData_->getRawFeatureData(featureIdx,*it));
  }

  datadefs::num_t leftFraction = 1.0*sampleIcs_left.size() / ( sampleIcs_left.size() + sampleIcs_right.size() );

  Node node;
  node.setSplitter(0,"foo",leftFraction,rawSplitValues_left,rawSplitValues_right);

  CPPUNIT_ASSERT( node.percolateData("1") == node.leftChild() );
  CPPUNIT_ASSERT( node.percolateData("2") == node.rightChild() );
  CPPUNIT_ASSERT( node.percolateData("3") == node.rightChild() );

  CPPUNIT_ASSERT( fabs( deltaImpurity - 0.001691905260604 ) < 1e-10 );

  

}

void treeDataTest::test_fullSplitterSweep() {

  Treedata* treeData = treeData_103by300NaN_;

  //size_t targetIdx = treeData->getFeatureIdx("N:output");
  size_t minSamples = 3;

  vector<string> checkList = utils::readListFromFile("test_fullSplitterSweep.txt",'\n');

  for ( size_t i = 0; i < checkList.size(); ++i ) {

    vector<string> checkLine = utils::split(checkList[i],'\t');

    CPPUNIT_ASSERT( checkLine.size() == 5 );

    size_t targetIdx = treeData->getFeatureIdx(checkLine[0]);
    size_t featureIdx = treeData->getFeatureIdx(checkLine[1]);

    CPPUNIT_ASSERT( featureIdx != targetIdx );
    CPPUNIT_ASSERT( featureIdx != treeData->end() );
    CPPUNIT_ASSERT( targetIdx != treeData->end() );

    num_t DI;
    size_t nIcs_left;
    size_t nIcs_right;

    vector<size_t> sampleIcs_left;
    vector<size_t> sampleIcs_right = utils::range( treeData->nSamples() );

    vector<num_t> tv = treeData->getFilteredFeatureData(targetIdx,sampleIcs_right);

    if ( treeData->isFeatureNumerical(featureIdx) ) {

      num_t splitValue;
      DI = treeData->numericalFeatureSplit(targetIdx,featureIdx,minSamples,sampleIcs_left,sampleIcs_right,splitValue);

    } else {
      
      set<num_t> splitValues_left, splitValues_right;
      DI = treeData->categoricalFeatureSplit(targetIdx,featureIdx,minSamples,sampleIcs_left,sampleIcs_right,splitValues_left,splitValues_right);
      
    }

    /*
      cout << treeData->getFeatureName(targetIdx) << " " << treeData->getFeatureName(featureIdx) << " " << utils::str2<num_t>(checkLine[2]) << " == " << DI 
      << " , " << utils::str2<size_t>(checkLine[3]) << " == " << sampleIcs_left.size() 
      << " , " << utils::str2<size_t>(checkLine[4]) << " == " << sampleIcs_right.size() << endl;
    */

    CPPUNIT_ASSERT( fabs( utils::str2<num_t>(checkLine[2]) - DI ) < 1e-5 );
    CPPUNIT_ASSERT( utils::str2<size_t>(checkLine[3]) == sampleIcs_left.size() );
    CPPUNIT_ASSERT( utils::str2<size_t>(checkLine[4]) == sampleIcs_right.size() );    

  }

}

void treeDataTest::test_end() {

  CPPUNIT_ASSERT( treeData_->getFeatureIdx("IDontExist") == treeData_->end() );

}

// Registers the fixture into the test 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION( treeDataTest ); 

#endif // TREEDATATEST_HPP
