#ifndef TREEDATA_NEWTEST_HPP
#define TREEDATA_NEWTEST_HPP

#include <cstdlib>

#include "newtest.hpp"

using namespace std;

void treedata_newtest_readAFM();
void treedata_newtest_readTransposedAFM();
void treedata_newtest_nRealSamples();

void treedata_newtest() {

  newtest( "Testing Treedata class with AFM data", &treedata_newtest_readAFM );
  newtest( "Testing Treedata class with transposed AFM data", &treedata_newtest_readTransposedAFM );
  newtest( "Testing proper counting of missing samples in Treedata", &treedata_newtest_nRealSamples );

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

void treedata_newtest_nRealSamples() {

  string fileName = "test/data/3by8_mixed_NA_matrix.afm";
  bool useContrasts = false;
  
  Treedata treeData(fileName,'\t',':',useContrasts);

  newassert(treeData.nFeatures() == 8);
  newassert(treeData.nRealSamples(0,1) == 1);
  newassert(treeData.nRealSamples(6,7) == 2);

}

#endif
