#ifndef TREEDATA_NEWTEST_HPP
#define TREEDATA_NEWTEST_HPP

#include <cstdlib>

#include "newtest.hpp"
#include "murmurhash3.hpp"
#include "distributions.hpp"

using namespace std;

void treedata_newtest_readAFM();
void treedata_newtest_readARFF();
void treedata_newtest_readTransposedAFM();
void treedata_newtest_nRealSamples();
void treedata_newtest_name2idxMap();
void treedata_newtest_numericalFeatureSplitsNumericalTarget();
void treedata_newtest_numericalFeatureSplitsCategoricalTarget();
void treedata_newtest_categoricalFeatureSplitsNumericalTarget();
void treedata_newtest_replaceFeatureData();
void treedata_newtest_end();
void treedata_newtest_hashFeature();
void treedata_newtest_bootstrapRealSamples();
void treedata_newtest_separateMissingSamples();

void treedata_newtest() {

  newtest( "Testing Treedata class with AFM data", &treedata_newtest_readAFM );
  newtest( "Testing Treedata clas with ARFF data", &treedata_newtest_readARFF );
  newtest( "Testing Treedata class with transposed AFM data", &treedata_newtest_readTransposedAFM );
  newtest( "Testing proper counting of missing samples in Treedata", &treedata_newtest_nRealSamples );
  newtest( "Testing mapping of feature names to indices and vice versa", &treedata_newtest_name2idxMap ); 
  newtest( "Testing splitting num->num", &treedata_newtest_numericalFeatureSplitsNumericalTarget );
  newtest( "Testing splitting num->cat", &treedata_newtest_numericalFeatureSplitsCategoricalTarget );
  newtest( "Testing splitting cat->num", &treedata_newtest_categoricalFeatureSplitsNumericalTarget );
  newtest( "Testing replacement of feature in Treedata", &treedata_newtest_replaceFeatureData );
  newtest( "Testing proper response for querying nonexistent feature" , &treedata_newtest_end );
  newtest( "Testing feature hashing", &treedata_newtest_hashFeature );
  newtest( "Testing bootstrapping", &treedata_newtest_bootstrapRealSamples );
  newtest( "Testing separation of missing samples", &treedata_newtest_separateMissingSamples );

}

void treedata_newtest_readAFM() {

  string fileName = "test/data/3by8_mixed_NA_matrix.afm";
  bool useContrasts = true;

  Reader reader(fileName,'\t');

  Treedata treeDataC(fileName,'\t',':',useContrasts);

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

  newassert( datadefs::isNAN(treeDataC.getFeatureData(0,0)) );
  newassert( fabs(treeDataC.getFeatureData(0,1) - 0.00) < datadefs::EPS );
  newassert( fabs(treeDataC.getFeatureData(0,2) - 0.000) < datadefs::EPS );
  newassert( treeDataC.getRawFeatureData(1,(size_t)0) == "foo" );
  newassert( datadefs::isNAN_STR(treeDataC.getRawFeatureData(1,(size_t)1)) );
  newassert( treeDataC.getRawFeatureData(1,(size_t)2) == "bar" );
  newassert( fabs(treeDataC.getFeatureData(2,0) - 2.2) < datadefs::EPS );
  newassert( fabs(treeDataC.getFeatureData(2,1) - 2.22) < datadefs::EPS );
  newassert( fabs(treeDataC.getFeatureData(2,2) - 2.222) < datadefs::EPS );
  newassert( fabs(treeDataC.getFeatureData(3,0) - 3.3) < datadefs::EPS );
  newassert( fabs(treeDataC.getFeatureData(3,1) - 3.33) < datadefs::EPS );
  newassert( fabs(treeDataC.getFeatureData(3,2) - 3.333) < datadefs::EPS );
  newassert( fabs(treeDataC.getFeatureData(4,0) - 4.4) < datadefs::EPS );
  newassert( fabs(treeDataC.getFeatureData(4,1) - 4.44) < datadefs::EPS );
  newassert( fabs(treeDataC.getFeatureData(4,2) - 4.444) < datadefs::EPS );
  newassert( fabs(treeDataC.getFeatureData(5,0) - 5.5) < datadefs::EPS );
  newassert( fabs(treeDataC.getFeatureData(5,1) - 5.55) < datadefs::EPS );
  newassert( fabs(treeDataC.getFeatureData(5,2) - 5.555) < datadefs::EPS );
  newassert( fabs(treeDataC.getFeatureData(6,0) - 6.6) < datadefs::EPS );
  newassert( datadefs::isNAN(treeDataC.getFeatureData(6,(size_t)1)) );
  newassert( fabs(treeDataC.getFeatureData(6,(size_t)2) - 6.666) < datadefs::EPS );

  useContrasts = false;
  
  Treedata treeData(fileName,'\t',':',useContrasts);

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

  newassert( datadefs::isNAN(treeData.getFeatureData(0,0)) );
  newassert( fabs(treeData.getFeatureData(0,1) - 0.00) < datadefs::EPS );
  newassert( fabs(treeData.getFeatureData(0,2) - 0.000) < datadefs::EPS );
  newassert( treeData.getRawFeatureData(1,(size_t)0) == "foo" );
  newassert( datadefs::isNAN_STR(treeData.getRawFeatureData(1,(size_t)1)) );
  newassert( treeData.getRawFeatureData(1,(size_t)2) == "bar" );
  newassert( fabs(treeData.getFeatureData(2,0) - 2.2) < datadefs::EPS );
  newassert( fabs(treeData.getFeatureData(2,1) - 2.22) < datadefs::EPS );
  newassert( fabs(treeData.getFeatureData(2,2) - 2.222) < datadefs::EPS );
  newassert( fabs(treeData.getFeatureData(3,0) - 3.3) < datadefs::EPS );
  newassert( fabs(treeData.getFeatureData(3,1) - 3.33) < datadefs::EPS );
  newassert( fabs(treeData.getFeatureData(3,2) - 3.333) < datadefs::EPS );
  newassert( fabs(treeData.getFeatureData(4,0) - 4.4) < datadefs::EPS );
  newassert( fabs(treeData.getFeatureData(4,1) - 4.44) < datadefs::EPS );
  newassert( fabs(treeData.getFeatureData(4,2) - 4.444) < datadefs::EPS );
  newassert( fabs(treeData.getFeatureData(5,0) - 5.5) < datadefs::EPS );
  newassert( fabs(treeData.getFeatureData(5,1) - 5.55) < datadefs::EPS );
  newassert( fabs(treeData.getFeatureData(5,2) - 5.555) < datadefs::EPS );
  newassert( fabs(treeData.getFeatureData(6,0) - 6.6) < datadefs::EPS );
  newassert( datadefs::isNAN(treeData.getFeatureData(6,1)) );
  newassert( fabs(treeData.getFeatureData(6,2) - 6.666) < datadefs::EPS );


}

void treedata_newtest_readTransposedAFM() {

  string fileName = "test/data/3by8_mixed_NA_transposed_matrix.afm";
  bool useContrasts = true;

  Treedata treeDataC(fileName,'\t',':',useContrasts);

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


  newassert( datadefs::isNAN(treeDataC.getFeatureData(0,0)) );
  newassert( fabs(treeDataC.getFeatureData(0,1) - 0.00) < datadefs::EPS );
  newassert( fabs(treeDataC.getFeatureData(0,2) - 0.000) < datadefs::EPS );
  newassert( treeDataC.getRawFeatureData(1,(size_t)0) == "foo" );
  newassert( datadefs::isNAN_STR(treeDataC.getRawFeatureData(1,(size_t)1)) );
  newassert( treeDataC.getRawFeatureData(1,(size_t)2) == "bar" );
  newassert( fabs(treeDataC.getFeatureData(2,0) - 2.2) < datadefs::EPS );
  newassert( fabs(treeDataC.getFeatureData(2,1) - 2.22) < datadefs::EPS );
  newassert( fabs(treeDataC.getFeatureData(2,2) - 2.222) < datadefs::EPS );
  newassert( fabs(treeDataC.getFeatureData(3,0) - 3.3) < datadefs::EPS );
  newassert( fabs(treeDataC.getFeatureData(3,1) - 3.33) < datadefs::EPS );
  newassert( fabs(treeDataC.getFeatureData(3,2) - 3.333) < datadefs::EPS );
  newassert( fabs(treeDataC.getFeatureData(4,0) - 4.4) < datadefs::EPS );
  newassert( fabs(treeDataC.getFeatureData(4,1) - 4.44) < datadefs::EPS );
  newassert( fabs(treeDataC.getFeatureData(4,2) - 4.444) < datadefs::EPS );
  newassert( fabs(treeDataC.getFeatureData(5,0) - 5.5) < datadefs::EPS );
  newassert( fabs(treeDataC.getFeatureData(5,1) - 5.55) < datadefs::EPS );
  newassert( fabs(treeDataC.getFeatureData(5,2) - 5.555) < datadefs::EPS );
  newassert( fabs(treeDataC.getFeatureData(6,0) - 6.6) < datadefs::EPS );
  newassert( datadefs::isNAN(treeDataC.getFeatureData(6,1)) );
  newassert( fabs(treeDataC.getFeatureData(6,2) - 6.666) < datadefs::EPS );


  useContrasts = false;

  Treedata treeData(fileName,'\t',':',useContrasts);

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

  newassert( datadefs::isNAN(treeData.getFeatureData(0,0)) );
  newassert( fabs(treeData.getFeatureData(0,1) - 0.00) < datadefs::EPS );
  newassert( fabs(treeData.getFeatureData(0,2) - 0.000) < datadefs::EPS );
  newassert( treeData.getRawFeatureData(1,(size_t)0) == "foo" );
  newassert( datadefs::isNAN_STR(treeData.getRawFeatureData(1,(size_t)1)) );
  newassert( treeData.getRawFeatureData(1,(size_t)2) == "bar" );
  newassert( fabs(treeData.getFeatureData(2,0) - 2.2) < datadefs::EPS );
  newassert( fabs(treeData.getFeatureData(2,1) - 2.22) < datadefs::EPS );
  newassert( fabs(treeData.getFeatureData(2,2) - 2.222) < datadefs::EPS );
  newassert( fabs(treeData.getFeatureData(3,0) - 3.3) < datadefs::EPS );
  newassert( fabs(treeData.getFeatureData(3,1) - 3.33) < datadefs::EPS );
  newassert( fabs(treeData.getFeatureData(3,2) - 3.333) < datadefs::EPS );
  newassert( fabs(treeData.getFeatureData(4,0) - 4.4) < datadefs::EPS );
  newassert( fabs(treeData.getFeatureData(4,1) - 4.44) < datadefs::EPS );
  newassert( fabs(treeData.getFeatureData(4,2) - 4.444) < datadefs::EPS );
  newassert( fabs(treeData.getFeatureData(5,0) - 5.5) < datadefs::EPS );
  newassert( fabs(treeData.getFeatureData(5,1) - 5.55) < datadefs::EPS );
  newassert( fabs(treeData.getFeatureData(5,2) - 5.555) < datadefs::EPS );
  newassert( fabs(treeData.getFeatureData(6,0) - 6.6) < datadefs::EPS );
  newassert( datadefs::isNAN(treeData.getFeatureData(6,1)) );
  newassert( fabs(treeData.getFeatureData(6,2) - 6.666) < datadefs::EPS );

}

void treedata_newtest_readARFF() {

  Treedata treeData("test_5by10_numeric_matrix.arff",'\t',':');

  newassert( treeData.nFeatures() == 5 );
  newassert( treeData.nSamples() == 10 );

  newassert( treeData.nRealSamples(0) == 10 );
  newassert( treeData.nRealSamples(1) == 9 );
  newassert( treeData.nRealSamples(2) == 9 );
  newassert( treeData.nRealSamples(3) == 10 );
  newassert( treeData.nRealSamples(4) == 10 );

  newassert( treeData.nRealSamples(0,1) == 9 );
  newassert( treeData.nRealSamples(1,2) == 8 );

  newassert( treeData.feature(0)->name() == "x1" );
  newassert( treeData.feature(1)->name() == "x2" );
  newassert( treeData.feature(2)->name() == "x3" );
  newassert( treeData.feature(3)->name() == "x4" );
  newassert( treeData.feature(4)->name() == "y"  );

  newassert( fabs( treeData.feature(0)->data[0] - 0.8147 ) < datadefs::EPS );
  newassert( fabs( treeData.feature(0)->data[1] - 0.9058 ) < datadefs::EPS );
  newassert( fabs( treeData.feature(0)->data[2] - 0.1270 ) < datadefs::EPS );
  newassert( fabs( treeData.feature(0)->data[3] - 0.9134 ) < datadefs::EPS );
  newassert( fabs( treeData.feature(0)->data[4] - 0.6324 ) < datadefs::EPS );
  newassert( fabs( treeData.feature(0)->data[5] - 0.0975 ) < datadefs::EPS );
  newassert( fabs( treeData.feature(0)->data[6] - 0.2785 ) < datadefs::EPS );
  newassert( fabs( treeData.feature(0)->data[7] - 0.5469 ) < datadefs::EPS );
  newassert( fabs( treeData.feature(0)->data[8] - 0.9575 ) < datadefs::EPS );
  newassert( fabs( treeData.feature(0)->data[9] - 0.9649 ) < datadefs::EPS );

  newassert( fabs( treeData.feature(1)->data[0] - 1.0000 ) < datadefs::EPS );
  newassert( fabs( treeData.feature(1)->data[1] - 2.0000 ) < datadefs::EPS );
  newassert( fabs( treeData.feature(1)->data[2] - 3.0000 ) < datadefs::EPS );
  newassert( fabs( treeData.feature(1)->data[3] - 4.0000 ) < datadefs::EPS );
  newassert( fabs( treeData.feature(1)->data[4] - 5.0000 ) < datadefs::EPS );
  newassert( fabs( treeData.feature(1)->data[5] - 6.0000 ) < datadefs::EPS );
  newassert( fabs( treeData.feature(1)->data[6] - 7.0000 ) < datadefs::EPS );
  newassert( datadefs::isNAN(treeData.feature(1)->data[7]) );
  newassert( fabs( treeData.feature(1)->data[8] - 9.0000 ) < datadefs::EPS );
  newassert( fabs( treeData.feature(1)->data[9] - 10.0000 ) < datadefs::EPS );

  newassert( fabs( treeData.feature(2)->data[0] - 0.0596 ) < datadefs::EPS );
  newassert( fabs( treeData.feature(2)->data[1] - 0.6820 ) < datadefs::EPS );
  newassert( fabs( treeData.feature(2)->data[2] - 0.0424 ) < datadefs::EPS );
  newassert( fabs( treeData.feature(2)->data[3] - 0.0714 ) < datadefs::EPS );
  newassert( datadefs::isNAN(treeData.feature(2)->data[4]) );
  newassert( fabs( treeData.feature(2)->data[5] - 0.0967 ) < datadefs::EPS );
  newassert( fabs( treeData.feature(2)->data[6] - 0.8181 ) < datadefs::EPS );
  newassert( fabs( treeData.feature(2)->data[7] - 0.8175 ) < datadefs::EPS );
  newassert( fabs( treeData.feature(2)->data[8] - 0.7224 ) < datadefs::EPS );
  newassert( fabs( treeData.feature(2)->data[9] - 0.1499 ) < datadefs::EPS );

  newassert( fabs( treeData.feature(3)->data[0] - 0.9160 ) < datadefs::EPS );
  newassert( fabs( treeData.feature(3)->data[1] - 0.0012 ) < datadefs::EPS );
  newassert( fabs( treeData.feature(3)->data[2] - 0.4624 ) < datadefs::EPS );
  newassert( fabs( treeData.feature(3)->data[3] - 0.4243 ) < datadefs::EPS );
  newassert( fabs( treeData.feature(3)->data[4] - 0.4609 ) < datadefs::EPS );
  newassert( fabs( treeData.feature(3)->data[5] - 0.7702 ) < datadefs::EPS );
  newassert( fabs( treeData.feature(3)->data[6] - 0.3225 ) < datadefs::EPS );
  newassert( fabs( treeData.feature(3)->data[7] - 0.7847 ) < datadefs::EPS );
  newassert( fabs( treeData.feature(3)->data[8] - 0.4714 ) < datadefs::EPS );
  newassert( fabs( treeData.feature(3)->data[9] - 0.0358 ) < datadefs::EPS );

  newassert( fabs( treeData.feature(4)->data[0] - 6.0000 ) < datadefs::EPS );
  newassert( fabs( treeData.feature(4)->data[1] - 14.0000 ) < datadefs::EPS );
  newassert( fabs( treeData.feature(4)->data[2] - 24.0000 ) < datadefs::EPS );
  newassert( fabs( treeData.feature(4)->data[3] - 36.0000 ) < datadefs::EPS );
  newassert( fabs( treeData.feature(4)->data[4] - 50.0000 ) < datadefs::EPS );
  newassert( fabs( treeData.feature(4)->data[5] - 66.0000 ) < datadefs::EPS );
  newassert( fabs( treeData.feature(4)->data[6] - 84.0000 ) < datadefs::EPS );
  newassert( fabs( treeData.feature(4)->data[7] - 104.0000 ) < datadefs::EPS );
  newassert( fabs( treeData.feature(4)->data[8] - 126.0000 ) < datadefs::EPS );
  newassert( fabs( treeData.feature(4)->data[9] - 150.0000 ) < datadefs::EPS );

}


void treedata_newtest_nRealSamples() {

  string fileName = "test/data/3by8_mixed_NA_matrix.afm";
  bool useContrasts = false;
  
  Treedata treeData(fileName,'\t',':',useContrasts);

  newassert(treeData.nFeatures() == 8);
  newassert(treeData.nRealSamples(0,1) == 1);
  newassert(treeData.nRealSamples(6,7) == 2);

}

void treedata_newtest_numericalFeatureSplitsNumericalTarget() {

  Treedata treeData("test_103by300_mixed_matrix.afm",'\t',':',true);

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

    newassert( fabs( deltaImpurity - 1.289680394982406 ) < 1e-10 );
    newassert( fabs( splitValue - 4.387 ) < 1e-10 );

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

  Treedata treeData("test_103by300_mixed_matrix.afm",'\t',':',true);

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

    newassert( fabs( deltaImpurity - 0.012389077212806 ) < 1e-10 );
    newassert( fabs( splitValue - 9.827 ) < 1e-10 );

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

  Treedata treeData("test_103by300_mixed_matrix.afm",'\t',':',true);

  vector<size_t> sampleIcs_left(0);
  vector<size_t> sampleIcs_right = utils::range(300);
  vector<size_t> sampleIcs_missing(0);

  unordered_set<num_t> splitValues_left,splitValues_right;
  datadefs::num_t deltaImpurity;

  size_t featureIdx = 1;
  size_t targetIdx = 0;
  size_t minSamples = 1;

  deltaImpurity = treeData.categoricalFeatureSplit(targetIdx,
						   featureIdx,
						   {1,2},
						   minSamples,
						   sampleIcs_left,
						   sampleIcs_right,
						   splitValues_left);
  

  //newassert( fabs( deltaImpurity - 1.102087375288799 ) < 1e-10 );

}

void treedata_newtest_end() {

  Treedata treeData("test_103by300_mixed_matrix.afm",'\t',':',true);

  newassert( treeData.getFeatureIdx("IDontExist") == treeData.end() );

}

void treedata_newtest_replaceFeatureData() {

  Treedata treeData("test_103by300_mixed_matrix.afm",'\t',':',true);

  treeData.replaceFeatureData(0,vector<num_t>(treeData.nSamples(),0.0));

  for ( size_t i = 0; i < treeData.nSamples(); ++i ) {
    newassert( treeData.getFeatureData(0,i) == 0.0 );
  }

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

  Treedata treeData("test_103by300_mixed_matrix.afm",'\t',':',true);

  distributions::Random random;

  bool withReplacement = true;
  num_t sampleSize = 1.0;
  size_t featureIdx = 0;
  vector<size_t> ics,oobIcs;

  num_t oobFraction = 0.0;

  for ( size_t i = 0; i < 1000; ++i ) {

    treeData.bootstrapFromRealSamples(&random,withReplacement,sampleSize,featureIdx,ics,oobIcs);

    oobFraction += 1.0 * oobIcs.size();

    newassert( ics.size() == treeData.nRealSamples(featureIdx) );
    newassert( !datadefs::containsNAN(treeData.getFeatureData(featureIdx,ics)) );
    newassert( !datadefs::containsNAN(treeData.getFeatureData(featureIdx,oobIcs)) );

  }

  oobFraction /= 1000 * treeData.nRealSamples(featureIdx);

  newassert( fabs( oobFraction - 0.36 ) < 0.05 );

}


void treedata_newtest_name2idxMap() {

  string fileName = "test_6by10_mixed_matrix.tsv";

  Treedata treeData(fileName,'\t',':',true);

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

  Treedata treeData(fileName,'\t',':');

  vector<size_t> sampleIcs = utils::range(10);
  vector<size_t> missingIcs;

  treeData.separateMissingSamples(0,sampleIcs,missingIcs);

  vector<num_t> filteredData = treeData.getFeatureData(0,sampleIcs);

  newassert( filteredData.size() == 8 );
  newassert( fabs( filteredData[0] - 8.5 ) < datadefs::EPS );
  newassert( fabs( filteredData[1] - 3.4 ) < datadefs::EPS );
  newassert( fabs( filteredData[2] - 7.2 ) < datadefs::EPS );
  newassert( fabs( filteredData[3] - 5 )   < datadefs::EPS );
  newassert( fabs( filteredData[4] - 6 )   < datadefs::EPS );
  newassert( fabs( filteredData[5] - 7 )   < datadefs::EPS );
  newassert( fabs( filteredData[6] - 11 )  < datadefs::EPS );
  newassert( fabs( filteredData[7] - 9 )   < datadefs::EPS );

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

  filteredData = treeData.getFeatureData(1,sampleIcs);

  newassert( filteredData.size() == 5 );
  newassert( fabs( filteredData[0] - 3 ) < datadefs::EPS );
  newassert( fabs( filteredData[1] - 4 ) < datadefs::EPS );
  newassert( fabs( filteredData[2] - 5 ) < datadefs::EPS );
  newassert( fabs( filteredData[3] - 6 ) < datadefs::EPS );
  newassert( fabs( filteredData[4] - 9 ) < datadefs::EPS );

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

  filteredData = treeData.getFeatureData(0,sampleIcs);

  newassert( filteredData.size() == 2 );
  newassert( fabs( filteredData[0] - 6 ) < datadefs::EPS );
  newassert( fabs( filteredData[1] - 6 ) < datadefs::EPS );

  newassert( sampleIcs.size() == 2 );
  newassert( sampleIcs[0] == 5 );
  newassert( sampleIcs[1] == 5 );

  treeData.separateMissingSamples(1,sampleIcs,missingIcs);

  filteredData = treeData.getFeatureData(1,sampleIcs);

  newassert( filteredData.size() == 0 );
  newassert( sampleIcs.size() == 0 );

}


#endif
