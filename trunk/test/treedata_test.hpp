#ifndef TREEDATATEST_HPP
#define TREEDATATEST_HPP

#include <cppunit/extensions/HelperMacros.h>
#include "treedata.hpp"
#include "datadefs.hpp"
#include "errno.hpp"

class treeDataTest : public CppUnit::TestFixture {
  CPPUNIT_TEST_SUITE( treeDataTest );
  //CPPUNIT_TEST( test_permuteContrasts );
  CPPUNIT_TEST( test_name2idxMap );
  CPPUNIT_TEST( test_getFeatureData );
  CPPUNIT_TEST( test_getRandomUnif );
  CPPUNIT_TEST( test_getFilteredFeatureData );
  CPPUNIT_TEST( test_getFilteredAndSortedFeatureDataPair3 );
  CPPUNIT_TEST( test_parseARFF );
  CPPUNIT_TEST( test_parseARFF_extended ); 
  CPPUNIT_TEST( test_keepFeatures );
  CPPUNIT_TEST( test_removeFeatures );
  CPPUNIT_TEST_SUITE_END();
  
public:
  void setUp();
  void tearDown();
  void test_permuteContrasts();
  void test_name2idxMap();
  void test_getFeatureData();
  void test_getRandomUnif();
  void test_getFilteredFeatureData();
  void test_getFilteredAndSortedFeatureDataPair3();
  void test_parseARFF();
  void test_parseARFF_extended();
  void test_keepFeatures();
  void test_removeFeatures();
};

void treeDataTest::setUp() {}
void treeDataTest::tearDown() {}
 

void treeDataTest::test_permuteContrasts() {
  
  string fileName = "test_6by10_mixed_matrix.tsv";
  
  Treedata::Treedata treeData(fileName,'\t',':');
  
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

  Treedata treedata("test_2by8_numerical_matrix.tsv",'\t',':');

  vector<size_t> sampleIcs = utils::range(8);

  vector<num_t> v1,v2;

  treedata.getFilteredFeatureDataPair(0,1,sampleIcs,v1,v2);

  

}


void treeDataTest::test_getFilteredFeatureData() {
  
  string fileName = "test_6by10_mixed_matrix.tsv";
  
  Treedata::Treedata treeData(fileName,'\t',':');


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

  Treedata::Treedata treeData(fileName,'\t',':');

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
  
  Treedata treeData("test_5by10_numeric_matrix.arff",'\t',':');

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

void treeDataTest::test_getRandomUnif() {

  Treedata treedata("test_6by10_mixed_matrix.tsv",'\t',':');

  for ( size_t i = 0; i < 100000; ++i ) {
    datadefs::num_t r = treedata.getRandomUnif();
    CPPUNIT_ASSERT( 0.0 <= r && r < 1.0 );
  }

}

void treeDataTest::test_parseARFF_extended() {

  Treedata treeData("test_12by21_categorical_matrix.arff",'\t',':');

}

void treeDataTest::test_keepFeatures() {
  
  Treedata treedata("test_6by10_mixed_matrix.tsv",'\t',':');

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

  Treedata treedata("test_6by10_mixed_matrix.tsv",'\t',':');

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


// Registers the fixture into the test 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION( treeDataTest ); 

#endif // TREEDATATEST_HPP
