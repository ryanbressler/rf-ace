#ifndef UTILS_NEWTEST_HPP
#define UTILS_NEWTEST_HPP

#include <cstdlib>
#include <vector>
#include <map>
#include <unordered_map>
#include <cmath>

#include "datadefs.hpp"
#include "newtest.hpp"
#include "utils.hpp"
#include "math.hpp"

using namespace std;
using datadefs::num_t;


void utils_newtest_categoricalFeatureSplitsNumericalTarget();
void utils_newtest_categoricalFeatureSplitsCategoricalTarget();
void utils_newtest_parse();
void utils_newtest_str2();
void utils_newtest_write();
void utils_newtest_range();
void utils_newtest_trim();
void utils_newtest_chomp();
void utils_newtest_split();
void utils_newtest_permute();
void utils_newtest_splitRange();
void utils_newtest_strv2catv();
void utils_newtest_strv2numv();
void utils_newtest_sortDataAndMakeRef();
void utils_newtest_sortFromRef();
void utils_newtest_text2tokens();



void utils_newtest() {

  newtest( "categoricalFeatureSplitsNumericalTarget(x)", &utils_newtest_categoricalFeatureSplitsNumericalTarget);
  newtest( "categoricalFeatureSplitsCategoricalTarget(x)", &utils_newtest_categoricalFeatureSplitsCategoricalTarget);
  newtest( "parse(x)", &utils_newtest_parse );
  newtest( "str2(x)", &utils_newtest_str2 );
  newtest( "write(x)", &utils_newtest_write );
  newtest( "range(x)", &utils_newtest_range );
  newtest( "trim(x)", &utils_newtest_trim );
  newtest( "chomp(x)", &utils_newtest_chomp );
  newtest( "split(x)", &utils_newtest_split );
  newtest( "permute(x)", &utils_newtest_permute );
  newtest( "splitRange(x)", &utils_newtest_splitRange );
  newtest( "strv2catv(x)", &utils_newtest_strv2catv );
  newtest( "strv2numv(x)", &utils_newtest_strv2numv );
  newtest( "sortAndMakeRef(x)", &utils_newtest_sortDataAndMakeRef );
  newtest( "sortFromRef(x)", &utils_newtest_sortFromRef );
  newtest( "text2tokens(x)", &utils_newtest_text2tokens );

}

void utils_newtest_categoricalFeatureSplitsNumericalTarget() {

  vector<num_t> fv = {1,1,1,2,2,2,3,3,3,4,4,4};
  vector<num_t> tv = {1,1,1,2,3,4,5,6,7,8,9,10};

  unordered_map<num_t,vector<size_t> > fmap_left,fmap_right;

  num_t DI = utils::categoricalFeatureSplitsNumericalTarget2(tv,fv,1,{1,2,3,4},fmap_left,fmap_right);
  
  num_t DI_ref = math::deltaImpurity_regr(math::mean(tv),12,math::mean({1,1,1,2,3,4}),6,math::mean({5,6,7,8,9,10}),6);
  
  newassert( fabs( DI - DI_ref ) < 1e-5 );

  fv = {1,1,1,1,1,1,1,1,1,1,1,1};

  DI = utils::categoricalFeatureSplitsNumericalTarget2(tv,fv,1,{1},fmap_left,fmap_right);

  DI_ref = 0;

  newassert( fabs( DI - DI_ref ) < 1e-3 );  
  
}

void utils_newtest_categoricalFeatureSplitsCategoricalTarget() {
  
  vector<num_t> fv = {1,1,1,2,2,2,3,3,3,4,4,4};
  vector<num_t> tv = {1,1,1,2,3,4,5,6,7,8,9,10};

  unordered_map<num_t,vector<size_t> > fmap_left,fmap_right;

  num_t DI = utils::categoricalFeatureSplitsCategoricalTarget2(tv,fv,1,{1,2,3,4},fmap_left,fmap_right);

  unordered_map<num_t,size_t> freq_left,freq_right,freq_tot;
  size_t sf_left = 0;
  size_t sf_right = 0;
  size_t sf_tot = 0;

  for ( size_t i = 0; i < tv.size(); ++i ) {
    math::incrementSquaredFrequency(tv[i],freq_tot,sf_tot);
  }
  
  for ( size_t i = 0; i < 3; ++i ) {
    math::incrementSquaredFrequency(tv[i],freq_left,sf_left);
  }

  for ( size_t i = 3; i < 12; ++i ) {
    math::incrementSquaredFrequency(tv[i],freq_right,sf_right);
  }
  
  num_t DI_ref = math::deltaImpurity_class(sf_tot,12,sf_left,3,sf_right,9);
 
  newassert( fabs( DI - DI_ref ) < 1e-5 );

  fv = {1,1,1,1,1,1,1,1,1,1,1,1};

  DI = utils::categoricalFeatureSplitsNumericalTarget2(tv,fv,1,{1},fmap_left,fmap_right);

  DI_ref = 0;

  newassert( fabs( DI - DI_ref ) < 1e-5 );

}

void utils_newtest_parse() {

  string s1("KEY1=val1,KEY2='val2',key3='val3,which=continues\"here'");

  map<string,string> m1 = utils::parse(s1,',','=','\'');

  newassert( m1["KEY1"] == "val1" );
  newassert( m1["KEY2"] == "val2" );
  newassert( m1["key3"] == "val3,which=continues\"here");

  string s2("KEY1=\"\",KEY2=,KEY3=\"=\"");

  map<string,string> m2 = utils::parse(s2,',','=','"');

  newassert( m2["KEY1"] == "" );
  newassert( m2["KEY2"] == "" );
  newassert( m2["KEY3"] == "=" );
  
}

void utils_newtest_str2() {

  string a("0.0");
  string b("1.0");
  string c("-1.0");
  string d("-1.0e10");

  newassert(utils::str2<num_t>(a) == 0.0);
  newassert(utils::str2<num_t>(b) == 1.0);
  newassert(utils::str2<num_t>(c) == -1.0);
  newassert(utils::str2<num_t>(d) == -1.0e10);

}

void utils_newtest_write() {

  vector<string> foo;
  foo.push_back("a");
  foo.push_back("b");

  stringstream is;
  string result;

  utils::write(is,foo.begin(),foo.end(),',');
  is >> result;
  is.clear();

  newassert( result == "a,b" );

  utils::write(is,foo.begin(),foo.end(),'a');
  is >> result;
  is.clear();

  newassert( result == "aab" );

  foo.pop_back();

  utils::write(is,foo.begin(),foo.end(),',');
  is >> result;
  is.clear();

  newassert( result == "a" );

  foo.pop_back();

  // Nothing should be read into "is", which is why "result" stays unchanged
  utils::write(is,foo.begin(),foo.end(),',');
  is >> result;
  is.clear();

  newassert( result == "a" );

}

void utils_newtest_range() {

  size_t n = 50;

  vector<size_t> ics = utils::range(n);

  for (size_t i = 0; i < n; ++i) {
    newassert(ics[i] == i);
  }
  
}

void utils_newtest_trim() {
  
  string wh1 = " \t";
  string wh2 = " ";

  string str = " \t  a \tb\t ";

  newassert( utils::trim(str,wh1) == "a \tb" );
  newassert( utils::trim(str,wh2) == "\t  a \tb\t" );

  str = "";
  
  newassert( utils::trim(str,wh1) == "" );
  newassert( utils::trim(str,wh2) == "" );

  str = "\t";

  newassert( utils::trim(str,wh1) == "" );
  newassert( utils::trim(str,wh2) == "\t" );

  str = " ";

  newassert( utils::trim(str,wh1) == "" );
  newassert( utils::trim(str,wh2) == "" );

}

void utils_newtest_chomp() {
  
  string str = "\tb\t\r";

  newassert( utils::chomp(str) == "\tb\t" );

  str = "\r  \t \r";

  newassert( utils::chomp(str) == "" );

  str = "";

  newassert( utils::chomp(str) == "" );

  str = "\r\n";

  newassert( utils::chomp(str) == "" );

  str = "a\n\r\n";

  newassert( utils::chomp(str) == "a" );

  str = "a\n";

  newassert( utils::chomp(str) == "a" );

}

void utils_newtest_split() {

  string str = " ab, c  , def,gh,,i j, ";

  vector<string> vec = utils::split(str,','," ");

  newassert( vec.size() == 7 );

  newassert( vec[0] == "ab" );
  newassert( vec[1] == "c" );
  newassert( vec[2] == "def" );
  newassert( vec[3] == "gh" );
  newassert( vec[4] == "" );
  newassert( vec[5] == "i j" );
  newassert( vec[6] == "" );

}

void utils_newtest_permute() {

  distributions::Random rand(0);

  num_t initData[] = {1.0,3.1,2.2,4.2,4.1,6.5,7.5,3,2};

  // Repeat the test 5 times
  for ( size_t i = 0; i < 5; ++i ) {
    
    vector<datadefs::num_t> data(initData,initData+8);
    vector<datadefs::num_t> dataOrig = data;
    
    utils::permute(data,&rand);
    
    bool anyChange = false;
    
    for ( size_t i = 0; i < data.size(); ++i ) {
      if ( data[i] != dataOrig[i] ) anyChange = true;
    }
    
    newassert( anyChange );
    
    sort(data.begin(),data.end());
    sort(dataOrig.begin(),dataOrig.end());
    
    for ( size_t i = 0; i < data.size(); ++i ) {
      newassert( data[i] == dataOrig[i] );
    }
  }

}

void utils_newtest_splitRange() {

  //vector<size_t> ics = utils::range(10);

  vector<vector<size_t> > icsSets = utils::splitRange(10,1);

  newassert( icsSets.size() == 1 );

  for ( size_t i = 0; i < icsSets.size(); ++i ) {
    for ( size_t j = 0; j < icsSets[i].size(); ++j ) {
      newassert( icsSets[i][j] == j );
    }
  }


}

void utils_newtest_strv2catv() {
  vector<string> strvec(51,"");
  vector<datadefs::num_t> catvec(51,0.0);
  map<string,datadefs::num_t> mapping;
  map<datadefs::num_t,string> backMapping;
  strvec[0] = "a";
  strvec[1] = "b";
  strvec[2] = "c";
  strvec[3] = "A";
  strvec[50] = "NaN";
  utils::strv2catv(strvec, catvec, mapping, backMapping);

  newassert(catvec[0] == 0.0);
  newassert(catvec[1] == 1.0);
  newassert(catvec[2] == 2.0);
  newassert(catvec[3] == 3.0);
  for (int i = 4; i < 50; ++i) {
    newassert(4.0);
  }
  newassert(datadefs::isNAN(catvec[50]));

  newassert( mapping["a"] == 0.0 );
  newassert( mapping["b"] == 1.0 );
  newassert( mapping["c"] == 2.0 );
  newassert( mapping["A"] == 3.0 );
  newassert( mapping[""]  == 4.0 );

  newassert( backMapping[0.0] == "a" );
  newassert( backMapping[1.0] == "b" );
  newassert( backMapping[2.0] == "c" );
  newassert( backMapping[3.0] == "A" );
  newassert( backMapping[4.0]  == "" );

}

void utils_newtest_strv2numv() {
  vector<string> strvec(51,"3.0");
  vector<datadefs::num_t> catvec(51,0.0);
  strvec[0] = "0.0";
  strvec[1] = "1.0";
  strvec[2] = "2.0";
  strvec[3] = "0.00";
  strvec[50] = "NaN";
  utils::strv2numv(strvec, catvec);

  newassert(catvec[0] == 0.0);
  newassert(catvec[1] == 1.0);
  newassert(catvec[2] == 2.0);
  newassert(catvec[3] == 0.0);
  for (int i = 4; i < 50; ++i) {
    newassert(catvec[i] == 3.0);
  }
  newassert(datadefs::isNAN(catvec[50]));
}

void utils_newtest_sortDataAndMakeRef() {
  vector<datadefs::num_t> data;
  vector<size_t> refIcs;
  for (int i = 0; i < 50; ++i) {
    data.push_back(static_cast<datadefs::num_t>(i));
  }

  for (int i = 49; i > -1; --i) { // Deliberately of length data.size() - 1
    refIcs.push_back(static_cast<size_t>(i));
  } // This allocation should be irrelevant. We keep it here deliberately to
    //  ensure the results are flattened in a safe manner; we expect a bounds
    //  checker to complain violently should that not be the case.

  utils::sortDataAndMakeRef(true, data, refIcs);
  for (int i = 0; i < 50; ++i) {
    newassert(data[i] == static_cast<datadefs::num_t>(i));
    newassert(refIcs[i] == static_cast<size_t>(i));
  }

  utils::sortDataAndMakeRef(false, data, refIcs);
  for (int i = 0; i < 50; ++i) {
    newassert(data[i] == static_cast<datadefs::num_t>(49-i));
    newassert(refIcs[i] == static_cast<size_t>(49-i));
  }
  newassert(data[49] == 0.0);
  newassert(refIcs[49] == 0);

  // Check for correct behavior with an empty data list and arbitrary refIcs
  data.clear();
  utils::sortDataAndMakeRef(true, data, refIcs);
  newassert(data.size() == 0);
  newassert(refIcs.size() == 0);

  utils::sortDataAndMakeRef(false, data, refIcs);
  newassert(data.size() == 0);
  newassert(refIcs.size() == 0);

  // NaNs are not checked as sorting targets, as their behavior is currently undefined
}

void utils_newtest_sortFromRef() {
  vector<int> data(50,0);
  vector<size_t> refIcs(50,0);
  for (int i = 0; i < 50; ++i) {
    data[i] = i;
    refIcs[i] = 49-i;
  }

  utils::sortFromRef<int>(data,refIcs);
  for (int i = 0; i < 50; ++i) {
    newassert(data[i] == 49-i);
  }
}

void utils_newtest_text2tokens() {

  //char text[] = "I want to, tokenizE  This!!.; it's so rad@";

  //vector<string> tokens = utils::text2tokens(text);

  string text = "I want to, tokenizE  This!!.; it's so rad";

  unordered_set<uint32_t> hashes = utils::hashText(text);
  unordered_set<uint32_t>::const_iterator it = hashes.begin();

  //newassert( hashes.size() == 9 );

  uint32_t h;

  MurmurHash3_x86_32("i",1,0,&h);
  newassert( hashes.find(h) != hashes.end() );

  MurmurHash3_x86_32("want",4,0,&h);
  newassert( hashes.find(h) != hashes.end() );

  MurmurHash3_x86_32("to",2,0,&h);
  newassert( hashes.find(h) != hashes.end() );

  MurmurHash3_x86_32("tokenize",8,0,&h);
  newassert( hashes.find(h) != hashes.end() );

  MurmurHash3_x86_32("this",4,0,&h);
  newassert( hashes.find(h) != hashes.end() );

  MurmurHash3_x86_32("it",2,0,&h);
  newassert( hashes.find(h) != hashes.end() );

  MurmurHash3_x86_32("s",1,0,&h);
  newassert( hashes.find(h) != hashes.end() );

  MurmurHash3_x86_32("so",2,0,&h);
  newassert( hashes.find(h) != hashes.end() );

  MurmurHash3_x86_32("rad",3,0,&h);
  newassert( hashes.find(h) != hashes.end() );


}

#endif
