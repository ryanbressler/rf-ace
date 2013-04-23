#ifndef TREEDATA_NEWTEST_HPP
#define TREEDATA_NEWTEST_HPP

#include <cstdlib>

#include "newtest.hpp"
#include "murmurhash3.hpp"
#include "distributions.hpp"
#include "densetreedata.hpp"

using namespace std;

void treedata_newtest_readAFM();
void treedata_newtest_readTransposedAFM();
void treedata_newtest_nRealSamples();
void treedata_newtest_name2idxMap();
void treedata_newtest_numericalFeatureSplitsNumericalTarget();
void treedata_newtest_numericalFeatureSplitsCategoricalTarget();
void treedata_newtest_categoricalFeatureSplitsNumericalTarget();
//void treedata_newtest_replaceFeatureData();
void treedata_newtest_end();
void treedata_newtest_hashFeature();
void treedata_newtest_bootstrapRealSamples();
void treedata_newtest_separateMissingSamples();

void treedata_newtest() {

  newtest( "readAFM(x)", &treedata_newtest_readAFM );
  newtest( "readTransposedAFM(x)", &treedata_newtest_readTransposedAFM );
  newtest( "nRealSamples(x)", &treedata_newtest_nRealSamples );
  newtest( "name2idxMap(x)", &treedata_newtest_name2idxMap ); 
  newtest( "numericalFeatureSplitsNumericalTarget(x)", &treedata_newtest_numericalFeatureSplitsNumericalTarget );
  newtest( "numericalFeatureSplitsCategoricalTarget(x)", &treedata_newtest_numericalFeatureSplitsCategoricalTarget );
  newtest( "categoricalFeatureSplitsNumericalTarget(x)", &treedata_newtest_categoricalFeatureSplitsNumericalTarget );
  //newtest( "replaceFeatureData(x)", &treedata_newtest_replaceFeatureData );
  newtest( "end(x)" , &treedata_newtest_end );
  newtest( "hashFeature(x)", &treedata_newtest_hashFeature );
  newtest( "bootstrapRealSamples(x)", &treedata_newtest_bootstrapRealSamples );
  newtest( "separateMissingSamples(x)", &treedata_newtest_separateMissingSamples );

}

void treedata_newtest_readAFM() {

  string fileName = "test/data/3by8_mixed_NA_matrix.afm";
  bool useContrasts = true;

  Reader reader(fileName,'\t');

  DenseTreeData treeDataC(fileName,'\t',':',useContrasts);

  newassert( treeDataC.nFeatures() == 8 );
  newassert( treeDataC.nSamples() == 3  );
  newassert( treeDataC.features_.size() == 16 );
  newassert( treeDataC.name2idx_.size() == 16 );

  reader.nextLine();
  reader.skipField();

  string featureName;
  reader >> featureName; newassert( treeDataC.feature(0)->name() == featureName ); newassert( featureName == "N:var0" ); newassert( treeDataC.getFeatureIdx("N:var0") == 0 );
  reader >> featureName; newassert( treeDataC.feature(1)->name() == featureName ); newassert( featureName == "C:var1" ); newassert( treeDataC.getFeatureIdx("C:var1") == 1 );
  reader >> featureName; newassert( treeDataC.feature(2)->name() == featureName ); newassert( featureName == "N:var2" ); newassert( treeDataC.getFeatureIdx("N:var2") == 2 );
  reader >> featureName; newassert( treeDataC.feature(3)->name() == featureName ); newassert( featureName == "N:var3" ); newassert( treeDataC.getFeatureIdx("N:var3") == 3 );
  reader >> featureName; newassert( treeDataC.feature(4)->name() == featureName ); newassert( featureName == "N:var4" ); newassert( treeDataC.getFeatureIdx("N:var4") == 4 );
  reader >> featureName; newassert( treeDataC.feature(5)->name() == featureName ); newassert( featureName == "N:var5" ); newassert( treeDataC.getFeatureIdx("N:var5") == 5 );
  reader >> featureName; newassert( treeDataC.feature(6)->name() == featureName ); newassert( featureName == "N:var6" ); newassert( treeDataC.getFeatureIdx("N:var6") == 6 );
  reader >> featureName; newassert( treeDataC.feature(7)->name() == featureName ); newassert( featureName == "T:var7" ); newassert( treeDataC.getFeatureIdx("T:var7") == 7 );

  newassert( treeDataC.feature(8)->name()   == "N:var0_CONTRAST" ); newassert( treeDataC.getFeatureIdx("N:var0_CONTRAST") ==  8 );
  newassert( treeDataC.feature(9)->name()   == "C:var1_CONTRAST" ); newassert( treeDataC.getFeatureIdx("C:var1_CONTRAST") ==  9 );
  newassert( treeDataC.feature(10)->name()  == "N:var2_CONTRAST" ); newassert( treeDataC.getFeatureIdx("N:var2_CONTRAST") == 10 );
  newassert( treeDataC.feature(11)->name()  == "N:var3_CONTRAST" ); newassert( treeDataC.getFeatureIdx("N:var3_CONTRAST") == 11 );
  newassert( treeDataC.feature(12)->name()  == "N:var4_CONTRAST" ); newassert( treeDataC.getFeatureIdx("N:var4_CONTRAST") == 12 );
  newassert( treeDataC.feature(13)->name()  == "N:var5_CONTRAST" ); newassert( treeDataC.getFeatureIdx("N:var5_CONTRAST") == 13 );
  newassert( treeDataC.feature(14)->name()  == "N:var6_CONTRAST" ); newassert( treeDataC.getFeatureIdx("N:var6_CONTRAST") == 14 );
  newassert( treeDataC.feature(15)->name()  == "T:var7_CONTRAST" ); newassert( treeDataC.getFeatureIdx("T:var7_CONTRAST") == 15 );

  newassert( treeDataC.feature(0)->isMissing(0) );
  newassert( fabs(treeDataC.feature(0)->getNumData(1) - 0.00) < 1e-5 );
  newassert( fabs(treeDataC.feature(0)->getNumData(2) - 0.000) < 1e-5 );
  newassert( treeDataC.feature(1)->getCatData(0) == "foo" );
  newassert( treeDataC.feature(1)->isMissing(1) );
  newassert( treeDataC.feature(1)->getCatData(2) == "bar" );
  newassert( fabs(treeDataC.feature(2)->getNumData(0) - 2.2) < 1e-5 );
  newassert( fabs(treeDataC.feature(2)->getNumData(1) - 2.22) < 1e-5 );
  newassert( fabs(treeDataC.feature(2)->getNumData(2) - 2.222) < 1e-5 );
  newassert( fabs(treeDataC.feature(3)->getNumData(0) - 3.3) < 1e-5 );
  newassert( fabs(treeDataC.feature(3)->getNumData(1) - 3.33) < 1e-5 );
  newassert( fabs(treeDataC.feature(3)->getNumData(2) - 3.333) < 1e-5 );
  newassert( fabs(treeDataC.feature(4)->getNumData(0) - 4.4) < 1e-5 );
  newassert( fabs(treeDataC.feature(4)->getNumData(1) - 4.44) < 1e-5 );
  newassert( fabs(treeDataC.feature(4)->getNumData(2) - 4.444) < 1e-5 );
  newassert( fabs(treeDataC.feature(5)->getNumData(0) - 5.5) < 1e-5 );
  newassert( fabs(treeDataC.feature(5)->getNumData(1) - 5.55) < 1e-5 );
  newassert( fabs(treeDataC.feature(5)->getNumData(2) - 5.555) < 1e-5 );
  newassert( fabs(treeDataC.feature(6)->getNumData(0) - 6.6) < 1e-5 );
  newassert( treeDataC.feature(6)->isMissing(1) );
  newassert( fabs(treeDataC.feature(6)->getNumData(2) - 6.666) < 1e-5 );

  useContrasts = false;
  
  DenseTreeData treeData(fileName,'\t',':',useContrasts);

  newassert( treeData.nFeatures() == 8 );
  newassert( treeData.nSamples() == 3 );
  newassert( treeData.features_.size() == 8 );
  newassert( treeData.name2idx_.size() == 8 );
  reader.rewind();
  reader.nextLine();
  reader.skipField();

  reader >> featureName; newassert( treeData.feature(0)->name() == featureName ); newassert( featureName == "N:var0" ); newassert( treeData.getFeatureIdx("N:var0") == 0 ); 
  reader >> featureName; newassert( treeData.feature(1)->name() == featureName ); newassert( featureName == "C:var1" ); newassert( treeData.getFeatureIdx("C:var1") == 1 );   
  reader >> featureName; newassert( treeData.feature(2)->name() == featureName ); newassert( featureName == "N:var2" ); newassert( treeData.getFeatureIdx("N:var2") == 2 );
  reader >> featureName; newassert( treeData.feature(3)->name() == featureName ); newassert( featureName == "N:var3" ); newassert( treeData.getFeatureIdx("N:var3") == 3 );
  reader >> featureName; newassert( treeData.feature(4)->name() == featureName ); newassert( featureName == "N:var4" ); newassert( treeData.getFeatureIdx("N:var4") == 4 );
  reader >> featureName; newassert( treeData.feature(5)->name() == featureName ); newassert( featureName == "N:var5" ); newassert( treeData.getFeatureIdx("N:var5") == 5 );
  reader >> featureName; newassert( treeData.feature(6)->name() == featureName ); newassert( featureName == "N:var6" ); newassert( treeData.getFeatureIdx("N:var6") == 6 );
  reader >> featureName; newassert( treeData.feature(7)->name() == featureName ); newassert( featureName == "T:var7" ); newassert( treeData.getFeatureIdx("T:var7") == 7 );

  newassert( treeData.feature(0)->isMissing(0) );
  newassert( fabs(treeData.feature(0)->getNumData(1) - 0.00) < 1e-5 );
  newassert( fabs(treeData.feature(0)->getNumData(2) - 0.000) < 1e-5 );
  newassert( treeData.feature(1)->getCatData(0) == "foo" );
  newassert( treeData.feature(1)->isMissing(1) );
  newassert( treeData.feature(1)->getCatData(2) == "bar" );
  newassert( fabs(treeData.feature(2)->getNumData(0) - 2.2) < 1e-5 );
  newassert( fabs(treeData.feature(2)->getNumData(1) - 2.22) < 1e-5 );
  newassert( fabs(treeData.feature(2)->getNumData(2) - 2.222) < 1e-5 );
  newassert( fabs(treeData.feature(3)->getNumData(0) - 3.3) < 1e-5 );
  newassert( fabs(treeData.feature(3)->getNumData(1) - 3.33) < 1e-5 );
  newassert( fabs(treeData.feature(3)->getNumData(2) - 3.333) < 1e-5 );
  newassert( fabs(treeData.feature(4)->getNumData(0) - 4.4) < 1e-5 );
  newassert( fabs(treeData.feature(4)->getNumData(1) - 4.44) < 1e-5 );
  newassert( fabs(treeData.feature(4)->getNumData(2) - 4.444) < 1e-5 );
  newassert( fabs(treeData.feature(5)->getNumData(0) - 5.5) < 1e-5 );
  newassert( fabs(treeData.feature(5)->getNumData(1) - 5.55) < 1e-5 );
  newassert( fabs(treeData.feature(5)->getNumData(2) - 5.555) < 1e-5 );
  newassert( fabs(treeData.feature(6)->getNumData(0) - 6.6) < 1e-5 );
  newassert( treeData.feature(6)->isMissing(1) );
  newassert( fabs(treeData.feature(6)->getNumData(2) - 6.666) < 1e-5 );

}

void treedata_newtest_readTransposedAFM() {

  string fileName = "test/data/3by8_mixed_NA_transposed_matrix.afm";
  bool useContrasts = true;

  DenseTreeData treeDataC(fileName,'\t',':',useContrasts);

  newassert( treeDataC.nFeatures() == 8 );
  newassert( treeDataC.nSamples() == 3  );
  newassert( treeDataC.features_.size() == 16 );
  newassert( treeDataC.name2idx_.size() == 16 );

  newassert( treeDataC.getFeatureIdx("N:var0") == 0 ); newassert( treeDataC.feature(0)->name() == "N:var0" );
  newassert( treeDataC.getFeatureIdx("C:var1") == 1 ); newassert( treeDataC.feature(1)->name() == "C:var1" );
  newassert( treeDataC.getFeatureIdx("N:var2") == 2 ); newassert( treeDataC.feature(2)->name() == "N:var2" );
  newassert( treeDataC.getFeatureIdx("N:var3") == 3 ); newassert( treeDataC.feature(3)->name() == "N:var3" );
  newassert( treeDataC.getFeatureIdx("N:var4") == 4 ); newassert( treeDataC.feature(4)->name() == "N:var4" );
  newassert( treeDataC.getFeatureIdx("N:var5") == 5 ); newassert( treeDataC.feature(5)->name() == "N:var5" );
  newassert( treeDataC.getFeatureIdx("N:var6") == 6 ); newassert( treeDataC.feature(6)->name() == "N:var6" );
  newassert( treeDataC.getFeatureIdx("T:var7") == 7 ); newassert( treeDataC.feature(7)->name() == "T:var7" );

  newassert( treeDataC.feature(8)->name()   == "N:var0_CONTRAST" ); newassert( treeDataC.getFeatureIdx("N:var0_CONTRAST") ==  8 );
  newassert( treeDataC.feature(9)->name()   == "C:var1_CONTRAST" ); newassert( treeDataC.getFeatureIdx("C:var1_CONTRAST") ==  9 );
  newassert( treeDataC.feature(10)->name()  == "N:var2_CONTRAST" ); newassert( treeDataC.getFeatureIdx("N:var2_CONTRAST") == 10 );
  newassert( treeDataC.feature(11)->name()  == "N:var3_CONTRAST" ); newassert( treeDataC.getFeatureIdx("N:var3_CONTRAST") == 11 );
  newassert( treeDataC.feature(12)->name()  == "N:var4_CONTRAST" ); newassert( treeDataC.getFeatureIdx("N:var4_CONTRAST") == 12 );
  newassert( treeDataC.feature(13)->name()  == "N:var5_CONTRAST" ); newassert( treeDataC.getFeatureIdx("N:var5_CONTRAST") == 13 );
  newassert( treeDataC.feature(14)->name()  == "N:var6_CONTRAST" ); newassert( treeDataC.getFeatureIdx("N:var6_CONTRAST") == 14 );
  newassert( treeDataC.feature(15)->name()  == "T:var7_CONTRAST" ); newassert( treeDataC.getFeatureIdx("T:var7_CONTRAST") == 15 );


  newassert( treeDataC.feature(0)->isMissing(0) );
  newassert( fabs(treeDataC.feature(0)->getNumData(1) - 0.00) < 1e-5 );
  newassert( fabs(treeDataC.feature(0)->getNumData(2) - 0.000) < 1e-5 );
  newassert( treeDataC.feature(1)->getCatData(0) == "foo" );
  newassert( treeDataC.feature(1)->isMissing(1) );
  newassert( treeDataC.feature(1)->getCatData(2) == "bar" );
  newassert( fabs(treeDataC.feature(2)->getNumData(0) - 2.2) < 1e-5 );
  newassert( fabs(treeDataC.feature(2)->getNumData(1) - 2.22) < 1e-5 );
  newassert( fabs(treeDataC.feature(2)->getNumData(2) - 2.222) < 1e-5 );
  newassert( fabs(treeDataC.feature(3)->getNumData(0) - 3.3) < 1e-5 );
  newassert( fabs(treeDataC.feature(3)->getNumData(1) - 3.33) < 1e-5 );
  newassert( fabs(treeDataC.feature(3)->getNumData(2) - 3.333) < 1e-5 );
  newassert( fabs(treeDataC.feature(4)->getNumData(0) - 4.4) < 1e-5 );
  newassert( fabs(treeDataC.feature(4)->getNumData(1) - 4.44) < 1e-5 );
  newassert( fabs(treeDataC.feature(4)->getNumData(2) - 4.444) < 1e-5 );
  newassert( fabs(treeDataC.feature(5)->getNumData(0) - 5.5) < 1e-5 );
  newassert( fabs(treeDataC.feature(5)->getNumData(1) - 5.55) < 1e-5 );
  newassert( fabs(treeDataC.feature(5)->getNumData(2) - 5.555) < 1e-5 );
  newassert( fabs(treeDataC.feature(6)->getNumData(0) - 6.6) < 1e-5 );
  newassert( treeDataC.feature(6)->isMissing(1) );
  newassert( fabs(treeDataC.feature(6)->getNumData(2) - 6.666) < 1e-5 );

  useContrasts = false;

  DenseTreeData treeData(fileName,'\t',':',useContrasts);

  newassert( treeData.nFeatures() == 8 );
  newassert( treeData.nSamples() == 3 );
  newassert( treeData.features_.size() == 8 );
  newassert( treeData.name2idx_.size() == 8 );

  newassert( treeData.getFeatureIdx("N:var0") == 0 ); newassert( treeData.feature(0)->name() == "N:var0" );
  newassert( treeData.getFeatureIdx("C:var1") == 1 ); newassert( treeData.feature(1)->name() == "C:var1" );
  newassert( treeData.getFeatureIdx("N:var2") == 2 ); newassert( treeData.feature(2)->name() == "N:var2" );
  newassert( treeData.getFeatureIdx("N:var3") == 3 ); newassert( treeData.feature(3)->name() == "N:var3" );
  newassert( treeData.getFeatureIdx("N:var4") == 4 ); newassert( treeData.feature(4)->name() == "N:var4" );
  newassert( treeData.getFeatureIdx("N:var5") == 5 ); newassert( treeData.feature(5)->name() == "N:var5" );
  newassert( treeData.getFeatureIdx("N:var6") == 6 ); newassert( treeData.feature(6)->name() == "N:var6" );
  newassert( treeData.getFeatureIdx("T:var7") == 7 ); newassert( treeData.feature(7)->name() == "T:var7" );

  newassert( treeData.getFeatureIdx("N:var0_CONTRAST") == treeData.end() );
  newassert( treeData.getFeatureIdx("C:var1_CONTRAST") == treeData.end() );
  newassert( treeData.getFeatureIdx("N:var2_CONTRAST") == treeData.end() );
  newassert( treeData.getFeatureIdx("N:var3_CONTRAST") == treeData.end() );
  newassert( treeData.getFeatureIdx("N:var4_CONTRAST") == treeData.end() );
  newassert( treeData.getFeatureIdx("N:var5_CONTRAST") == treeData.end() );
  newassert( treeData.getFeatureIdx("N:var6_CONTRAST") == treeData.end() );
  newassert( treeData.getFeatureIdx("T:var7_CONTRAST") == treeData.end() );

  newassert( treeData.feature(0)->isMissing(0) );
  newassert( fabs(treeData.feature(0)->getNumData(1) - 0.00) < 1e-5 );
  newassert( fabs(treeData.feature(0)->getNumData(2) - 0.000) < 1e-5 );
  newassert( treeData.feature(1)->getCatData(0) == "foo" );
  newassert( treeData.feature(1)->isMissing(1) );
  newassert( treeData.feature(1)->getCatData(2) == "bar" );
  newassert( fabs(treeData.feature(2)->getNumData(0) - 2.2) < 1e-5 );
  newassert( fabs(treeData.feature(2)->getNumData(1) - 2.22) < 1e-5 );
  newassert( fabs(treeData.feature(2)->getNumData(2) - 2.222) < 1e-5 );
  newassert( fabs(treeData.feature(3)->getNumData(0) - 3.3) < 1e-5 );
  newassert( fabs(treeData.feature(3)->getNumData(1) - 3.33) < 1e-5 );
  newassert( fabs(treeData.feature(3)->getNumData(2) - 3.333) < 1e-5 );
  newassert( fabs(treeData.feature(4)->getNumData(0) - 4.4) < 1e-5 );
  newassert( fabs(treeData.feature(4)->getNumData(1) - 4.44) < 1e-5 );
  newassert( fabs(treeData.feature(4)->getNumData(2) - 4.444) < 1e-5 );
  newassert( fabs(treeData.feature(5)->getNumData(0) - 5.5) < 1e-5 );
  newassert( fabs(treeData.feature(5)->getNumData(1) - 5.55) < 1e-5 );
  newassert( fabs(treeData.feature(5)->getNumData(2) - 5.555) < 1e-5 );
  newassert( fabs(treeData.feature(6)->getNumData(0) - 6.6) < 1e-5 );
  newassert( treeData.feature(6)->isMissing(1) );
  newassert( fabs(treeData.feature(6)->getNumData(2) - 6.666) < 1e-5 );

}

void treedata_newtest_nRealSamples() {

  string fileName = "test/data/3by8_mixed_NA_matrix.afm";
  bool useContrasts = false;
  
  DenseTreeData treeData(fileName,'\t',':',useContrasts);

  newassert(treeData.nFeatures() == 8);

}

void treedata_newtest_numericalFeatureSplitsNumericalTarget() {

  DenseTreeData treeData("test_103by300_mixed_matrix.afm",'\t',':',true);

  vector<size_t> sampleIcs_left(0);
  vector<size_t> sampleIcs_right = utils::range(300);
  vector<size_t> sampleIcs_missing(0);

  datadefs::num_t splitValue;
  datadefs::num_t deltaImpurity;

  size_t minSamples = 1;

  size_t targetIdx = 0; // numerical
  size_t featureIdx = 2; // numerical

  deltaImpurity = treeData.numericalFeatureSplit(targetIdx,
						 featureIdx,
						 minSamples,
						 sampleIcs_left,
						 sampleIcs_right,
						 splitValue);
  
  {
    set<size_t> leftIcs(sampleIcs_left.begin(),sampleIcs_left.end());
    set<size_t> rightIcs(sampleIcs_right.begin(),sampleIcs_right.end());

    newassert( fabs( deltaImpurity - 1.289680394982406 ) < 1e-5 );
    newassert( fabs( splitValue - 4.387 ) < 1e-5 );

    newassert( sampleIcs_left.size() == 127 );
    newassert( sampleIcs_right.size() == 173 );

    newassert( leftIcs.find(198) != leftIcs.end() );
    newassert( leftIcs.find(8)   != leftIcs.end() );
    newassert( leftIcs.find(102) != leftIcs.end() );
    newassert( leftIcs.find(4)   != leftIcs.end() );
    newassert( leftIcs.find(299) != leftIcs.end() );

    newassert( rightIcs.find(26) != rightIcs.end() );
    newassert( rightIcs.find(2)  != rightIcs.end() );
    newassert( rightIcs.find(10) != rightIcs.end() );
    newassert( rightIcs.find(81) != rightIcs.end() );
    newassert( rightIcs.find(33) != rightIcs.end() );
  }

  minSamples = 50;

}

void treedata_newtest_numericalFeatureSplitsCategoricalTarget() {

  DenseTreeData treeData("test_103by300_mixed_matrix.afm",'\t',':',true);

  size_t targetIdx = 1; // categorical
  size_t featureIdx = 2; // numerical

  vector<size_t> sampleIcs_left(0);
  vector<size_t> sampleIcs_right = utils::range(300);
  vector<size_t> sampleIcs_missing(0);

  datadefs::num_t splitValue;
  datadefs::num_t deltaImpurity;

  size_t minSamples = 1;

  deltaImpurity = treeData.numericalFeatureSplit(targetIdx,
						  featureIdx,
						  minSamples,
						  sampleIcs_left,
						  sampleIcs_right,
						  splitValue);
  
  {
    set<size_t> leftIcs(sampleIcs_left.begin(),sampleIcs_left.end());
    set<size_t> rightIcs(sampleIcs_right.begin(),sampleIcs_right.end());

    newassert( fabs( deltaImpurity - 0.012389077212806 ) < 1e-5 );
    newassert( fabs( splitValue - 9.827 ) < 1e-5 );

    newassert( sampleIcs_left.size() == 295 );
    newassert( sampleIcs_right.size() == 5 );

    newassert( leftIcs.find(261) != leftIcs.end() );
    newassert( leftIcs.find(185) != leftIcs.end() );
    newassert( leftIcs.find(3)   != leftIcs.end() );
    newassert( leftIcs.find(7)   != leftIcs.end() );
    newassert( leftIcs.find(256) != leftIcs.end() );

    newassert( rightIcs.find(69)  != rightIcs.end() );
    newassert( rightIcs.find(55)  != rightIcs.end() );
    newassert( rightIcs.find(100) != rightIcs.end() );
    newassert( rightIcs.find(127) != rightIcs.end() );
    newassert( rightIcs.find(91)  != rightIcs.end() );
  }
}

void treedata_newtest_categoricalFeatureSplitsNumericalTarget() {

  DenseTreeData treeData("test_103by300_mixed_matrix.afm",'\t',':',true);

  vector<size_t> sampleIcs_left(0);
  vector<size_t> sampleIcs_right = utils::range(300);
  vector<size_t> sampleIcs_missing(0);

  unordered_set<cat_t> splitValues_left,splitValues_right;
  
  size_t featureIdx = 1;
  size_t targetIdx = 0;
  size_t minSamples = 1;
  
  datadefs::num_t deltaImpurity = treeData.categoricalFeatureSplit(targetIdx,
								   featureIdx,
								   {"1","2"},
								   minSamples,
								   sampleIcs_left,
								   sampleIcs_right,
								   splitValues_left);
  

  newassert( fabs( deltaImpurity - 1.102087375288799 ) < 1e-5 );

}

void treedata_newtest_end() {

  DenseTreeData treeData("test_103by300_mixed_matrix.afm",'\t',':',true);

  newassert( treeData.getFeatureIdx("IDontExist") == treeData.end() );

}

void treedata_newtest_hashFeature() {

  vector<string> textData(3,"");

  textData[0] = "I am a random, text that is going to be hashed!";
  textData[1] = "I am a another random text";
  textData[2] = "This is something, completely different";

  bool doHash = true;

  Feature hashFeature(textData,"T:foo",doHash);

  uint32_t h;

  MurmurHash3_x86_32("i",1,0,&h);
  newassert( hashFeature.hasHash(0,h) );
  newassert( hashFeature.hasHash(1,h) );
  newassert( ! hashFeature.hasHash(2,h) );

  MurmurHash3_x86_32("am",2,0,&h);
  newassert( hashFeature.hasHash(0,h) );
  newassert( hashFeature.hasHash(1,h) );
  newassert( ! hashFeature.hasHash(2,h) );

  MurmurHash3_x86_32("random",6,0,&h);
  newassert(   hashFeature.hasHash(0,h) );
  newassert(   hashFeature.hasHash(1,h) );
  newassert( ! hashFeature.hasHash(2,h) );

  MurmurHash3_x86_32("text",4,0,&h);
  newassert(   hashFeature.hasHash(0,h) );

}

void treedata_newtest_bootstrapRealSamples() {

  DenseTreeData treeData("test_103by300_mixed_matrix.afm",'\t',':',true);

  distributions::Random random;

  bool withReplacement = true;
  num_t sampleSize = 1.0;
  size_t featureIdx = 0;
  vector<size_t> ics,oobIcs;

  num_t oobFraction = 0.0;

  for ( size_t i = 0; i < 1000; ++i ) {

    treeData.bootstrapFromRealSamples(&random,withReplacement,sampleSize,featureIdx,ics,oobIcs);

    oobFraction += 1.0 * oobIcs.size();

    newassert( ics.size() == treeData.feature(featureIdx)->nRealSamples() );
    newassert( !datadefs::containsNAN(treeData.feature(featureIdx)->getNumData(ics)) );
    newassert( !datadefs::containsNAN(treeData.feature(featureIdx)->getNumData(oobIcs)) );

  }

  oobFraction /= 1000 * treeData.feature(featureIdx)->nRealSamples();

  newassert( fabs( oobFraction - 0.36 ) < 0.05 );

}


void treedata_newtest_name2idxMap() {

  string fileName = "test_6by10_mixed_matrix.tsv";

  DenseTreeData treeData(fileName,'\t',':',true);

  newassert( treeData.name2idx_.size() == 2*treeData.nFeatures() );
  newassert( treeData.features_.size() == 2*treeData.nFeatures() );

  newassert( treeData.feature(0)->name() == "N:F1" );
  newassert( treeData.feature(1)->name() == "N:F2" );
  newassert( treeData.feature(2)->name() == "C:F3" );
  newassert( treeData.feature(3)->name() == "N:F4" );
  newassert( treeData.feature(4)->name() == "C:F5" );
  newassert( treeData.feature(5)->name() == "N:F6" );

  newassert( treeData.feature(6)->name() == "N:F1_CONTRAST" );
  newassert( treeData.feature(7)->name() == "N:F2_CONTRAST" );
  newassert( treeData.feature(8)->name() == "C:F3_CONTRAST" );
  newassert( treeData.feature(9)->name() == "N:F4_CONTRAST" );
  newassert( treeData.feature(10)->name() == "C:F5_CONTRAST" );
  newassert( treeData.feature(11)->name() == "N:F6_CONTRAST" );

  newassert( treeData.name2idx_["N:F1"] == 0 );
  newassert( treeData.name2idx_["N:F2"] == 1 );
  newassert( treeData.name2idx_["C:F3"] == 2 );
  newassert( treeData.name2idx_["N:F4"] == 3 );
  newassert( treeData.name2idx_["C:F5"] == 4 );
  newassert( treeData.name2idx_["N:F6"] == 5 );

  newassert( treeData.name2idx_["N:F1_CONTRAST"] == 6 );
  newassert( treeData.name2idx_["N:F2_CONTRAST"] == 7 );
  newassert( treeData.name2idx_["C:F3_CONTRAST"] == 8 );
  newassert( treeData.name2idx_["N:F4_CONTRAST"] == 9 );
  newassert( treeData.name2idx_["C:F5_CONTRAST"] == 10 );
  newassert( treeData.name2idx_["N:F6_CONTRAST"] == 11 );

}

void treedata_newtest_separateMissingSamples() {

  string fileName = "test_6by10_mixed_matrix.tsv";

  DenseTreeData treeData(fileName,'\t',':');

  vector<size_t> sampleIcs = utils::range(10);
  vector<size_t> missingIcs;

  treeData.separateMissingSamples(0,sampleIcs,missingIcs);

  vector<num_t> filteredData = treeData.feature(0)->getNumData(sampleIcs);

  newassert( filteredData.size() == 8 );
  newassert( fabs( filteredData[0] - 8.5 ) < 1e-5 );
  newassert( fabs( filteredData[1] - 3.4 ) < 1e-5 );
  newassert( fabs( filteredData[2] - 7.2 ) < 1e-5 );
  newassert( fabs( filteredData[3] - 5 )   < 1e-5 );
  newassert( fabs( filteredData[4] - 6 )   < 1e-5 );
  newassert( fabs( filteredData[5] - 7 )   < 1e-5 );
  newassert( fabs( filteredData[6] - 11 )  < 1e-5 );
  newassert( fabs( filteredData[7] - 9 )   < 1e-5 );

  newassert( sampleIcs.size() == 8 );

  newassert( sampleIcs[1] == 2 );
  newassert( sampleIcs[2] == 3 );
  newassert( sampleIcs[3] == 4 );
  newassert( sampleIcs[4] == 5 );
  newassert( sampleIcs[5] == 6 );
  newassert( sampleIcs[6] == 7 );
  newassert( sampleIcs[7] == 8 );

  //sampleIcs = utils::range(10);
  treeData.separateMissingSamples(1,sampleIcs,missingIcs);

  filteredData = treeData.feature(1)->getNumData(sampleIcs);

  newassert( filteredData.size() == 5 );
  newassert( fabs( filteredData[0] - 3 ) < 1e-5 );
  newassert( fabs( filteredData[1] - 4 ) < 1e-5 );
  newassert( fabs( filteredData[2] - 5 ) < 1e-5 );
  newassert( fabs( filteredData[3] - 6 ) < 1e-5 );
  newassert( fabs( filteredData[4] - 9 ) < 1e-5 );

  newassert( sampleIcs.size() == 5 );
  newassert( sampleIcs[0] == 1 );
  newassert( sampleIcs[1] == 2 );
  newassert( sampleIcs[2] == 3 );
  newassert( sampleIcs[3] == 4 );
  newassert( sampleIcs[4] == 7 );

  sampleIcs[0] = 0;
  sampleIcs[1] = 0;
  sampleIcs[2] = 0;
  sampleIcs[3] = 5;
  sampleIcs[4] = 5;

  treeData.separateMissingSamples(0,sampleIcs,missingIcs);

  filteredData = treeData.feature(0)->getNumData(sampleIcs);

  newassert( filteredData.size() == 2 );
  newassert( fabs( filteredData[0] - 6 ) < 1e-5 );
  newassert( fabs( filteredData[1] - 6 ) < 1e-5 );

  newassert( sampleIcs.size() == 2 );
  newassert( sampleIcs[0] == 5 );
  newassert( sampleIcs[1] == 5 );

  treeData.separateMissingSamples(1,sampleIcs,missingIcs);

  filteredData = treeData.feature(1)->getNumData(sampleIcs);

  newassert( filteredData.size() == 0 );
  newassert( sampleIcs.size() == 0 );
  
}


#endif
