#ifndef TREEDATA_NEWTEST_HPP
#define TREEDATA_NEWTEST_HPP

#include <cstdlib>

#include "newtest.hpp"

using namespace std;

void treedata_newtest_readAFM();

void treedata_newtest() {

  newtest( "Testing Treedata class with AFM data", &treedata_newtest_readAFM );

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

}

#endif
