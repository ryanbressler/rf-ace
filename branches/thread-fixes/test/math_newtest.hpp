#ifndef MATH_NEWTEST_HPP
#define MATH_NEWTEST_HPP

#include <cmath>
#include "math.hpp"
#include "datadefs.hpp"
#include "newtest.hpp"

using namespace std;

void math_newtest_erf();
void math_newtest_var();
void math_newtest_pearsonCorrelation();
void math_newtest_ttest();
void math_newtest_mean();
void math_newtest_mode();
void math_newtest_gamma();
void math_newtest_incrementDecrementSquaredFrequency();
void math_newtest_deltaImpurity_class();

void math_newtest() {

  newtest( "erf(x)", &math_newtest_erf );
  newtest( "var(x)", &math_newtest_var );
  newtest( "corr(x)", &math_newtest_pearsonCorrelation );
  newtest( "ttest(x,y)", &math_newtest_ttest );
  newtest( "mean(x)", &math_newtest_mean );
  newtest( "mode(x)", &math_newtest_mode );
  newtest( "gamma(x)", &math_newtest_gamma );
  newtest( "sqFreq(x)", &math_newtest_incrementDecrementSquaredFrequency );
  newtest( "deltaImpurity_class(x)", &math_newtest_deltaImpurity_class );

}

void math_newtest_erf() {

  datadefs::num_t x;
  
  x = 0.0;
  newassert(math::erf(x) == 0.0);
  
  x = 1.0;
  newassert(fabs(math::erf(x) - 0.84292558224862812465971728670410811901092529296875) < 1e-5 );
  
  x = datadefs::NUM_NAN;
  newassert(datadefs::isNAN(math::erf(x)));

}

void math_newtest_var() {

  vector<datadefs::num_t> data;
  datadefs::num_t mu, var;
  //size_t nRealValues = static_cast<size_t>(-1);
  for (int i = 0; i < 50; ++i) {
    data.push_back(static_cast<datadefs::num_t>(i));
  }

  var = math::var(data);

  mu = math::mean(data);

  newassert( fabs( var - 212.5 ) < 1e-3 );

  var = math::var(data,mu);

  newassert( fabs( var - 212.5 ) < 1e-3 );

  data.push_back(datadefs::NUM_NAN);

  newassert( datadefs::isNAN(math::var(data)) );

  data.resize(0);

  newassert( datadefs::isNAN(math::var(data)) );

  data.resize(1,0.0);

  newassert( datadefs::isNAN(math::var(data)) );

}

void math_newtest_ttest() {

  vector<datadefs::num_t> x,y;

  for(size_t i = 0; i < 10; ++i) {
    x.push_back(static_cast<datadefs::num_t>(i+6.0));
    y.push_back(static_cast<datadefs::num_t>(i+4.0));
  }

  // While increasing larger and larger samples, we expect
  // p-value to decrease gradually
  newassert( fabs( math::ttest(x,y) - 0.078466097923427 ) < 1e-5 );
  newassert( fabs( math::ttest(y,x) - 0.921533902076573 ) < 1e-5 );
  newassert( fabs( math::ttest(x,y,true) - 0.078466097923427 ) < 1e-5 );
  newassert( fabs( math::ttest(y,x,true) - 0.921533902076573 ) < 1e-5 );
  x.push_back(16.0);
  newassert( fabs( math::ttest(x,y) - 0.044077929403175 ) < 1e-5 );
  newassert( fabs( math::ttest(y,x) - 0.955922070596825 ) < 1e-5 );
  newassert( fabs( math::ttest(x,y,true) - 0.043411097004321 ) < 1e-5 );
  newassert( fabs( math::ttest(y,x,true) - 0.956588902995679 ) < 1e-5 );
  y.pop_back();
  newassert( fabs( math::ttest(x,y) - 0.021736196542511 ) < 1e-5 );
  newassert( fabs( math::ttest(y,x) - 0.978263803457489 ) < 1e-5 );
  newassert( fabs( math::ttest(x,y,true) - 0.019925334730836 ) < 1e-5 );
  newassert( fabs( math::ttest(y,x,true) - 0.980074665269164 ) < 1e-5 );
  x.push_back(17.0);
  newassert( fabs( math::ttest(x,y) - 0.01 ) < 0.01 );
  x.push_back(18.0);
  newassert( fabs( math::ttest(x,y) - 0.007 ) < 0.01 );
  x.push_back(19.0);
  newassert( fabs( math::ttest(x,y) - 0.005 ) < 0.01 );
  x.push_back(20.0);
  newassert( fabs( math::ttest(x,y) - 0.003 ) < 0.01 );
  y.pop_back();
  newassert( fabs( math::ttest(x,y) - 0.002 ) < 0.01 );

  datadefs::num_t p_old = 1.0;

  // While increasing values in x, we expect p-value to increase
  for ( size_t i = 0; i < 100; ++i ) {
    for ( size_t j = 0; j < x.size(); ++j ) {
      ++x[j];
    }

    datadefs::num_t p_new = math::ttest(x,y);
    //cout << "TEST " << i << ": " << p_new << " <= " << p_old << endl;
    newassert( p_new <= p_old );
    newassert( p_new < 0.002 );
    p_old = p_new;
  }

  // Taking the t-test backwards should yield a large p-value
  newassert( math::ttest(y,x) > 0.999 );

  p_old = 0.0;

  // While increasing values in y, we expect p-values to increase
  for ( size_t i = 0; i < 30; ++i ) {
    for ( size_t j = 0; j < y.size(); ++j ) {
      y[j] += 2;
    }

    datadefs::num_t p_new = math::ttest(x,y);
    //cout << "TEST " << i << ": " << p_new << " <= " << p_old << endl;
    newassert( p_new >= p_old );
    p_old = p_new;
  }

  // We should now have a large p-value in p_old
  //newassert( p_old > 0.999 );

  x.clear();
  x.resize(20,5.0);
  y.clear();
  y.resize(20,6.0);

  // p-value with constant x smaller than y should be 1.0
  newassert( fabs( math::ttest(x,y) - 1.0 ) < datadefs::EPS );

  y.clear();
  y.resize(20,4.0);

  // p-value with constant x larger than y should be EPS
  newassert( fabs( math::ttest(x,y) ) == datadefs::EPS );

  y.clear();
  y.resize(20,5.0);

  // p-value with equal constant data should be 0.5
  newassert( fabs( math::ttest(x,y) - 0.5 ) < datadefs::EPS );

  p_old = 0.5;

  // Start increasing x slightly to make sure p-values behave smoothly
  for ( size_t i = 0; i < 100; ++ i) {
    x[0] += 10*datadefs::EPS;

    datadefs::num_t p_new = math::ttest(x,y);

    newassert( p_new <= p_old );

    p_old = p_new;
  }

  x[0] = 5.0;

  p_old = 0.5;

  // Start decreasing x slightly to make sure p-values behave smoothly
  for ( size_t i = 0; i < 100; ++ i) {
    x[0] -= 10*datadefs::EPS;

    datadefs::num_t p_new = math::ttest(x,y);

    newassert( p_new >= p_old );

    p_old = p_new;
  }

  x.clear();
  y.clear();

  for( size_t i = 0; i < 10; ++i ) {
    x.push_back(static_cast<datadefs::num_t>(i));
    y.push_back(static_cast<datadefs::num_t>(i));
  }

  x[9] = static_cast<datadefs::num_t>(9.1);

  //cout << math::ttest(x,y) << endl;

  newassert( fabs( math::ttest(x,y) - 0.497 ) < 0.01 );

}

void math_newtest_mean() {

  vector<datadefs::num_t> x;

  newassert( datadefs::isNAN(math::mean(x)) );
  x.push_back(1.0);
  newassert( fabs( math::mean(x) - 1.0 ) < datadefs::EPS );
  x.push_back(2.0);
  newassert( fabs( math::mean(x) - 1.5 ) < datadefs::EPS );
  x.push_back(datadefs::NUM_NAN);
  newassert( datadefs::isNAN(math::mean(x)) );

}

void math_newtest_mode() {

  vector<datadefs::num_t> x(20,0);
  newassert( math::mode(x) == 0 );

  vector<size_t> y = utils::range(20);

  y[19] = 18;
  newassert( math::mode(y) == 18 );

  y[18] = 17;
  newassert( math::mode(y) == 17 );

}

void math_newtest_gamma() {

  vector<datadefs::num_t> data;
  size_t numClasses = 5;

  for (int i = 0; i < 50; ++i) {
    data.push_back(static_cast<datadefs::num_t>(i));
  }

  datadefs::num_t leafGamma = math::gamma(data, numClasses);
  newassert( fabs(leafGamma + 0.025) < 1e-5 );

}

void math_newtest_pearsonCorrelation() {
  vector<datadefs::num_t> x;
  vector<datadefs::num_t> y;
  vector<datadefs::num_t> y2;
  datadefs::num_t corr;

  for (int i = 0; i < 50; ++i) {
    x.push_back(static_cast<datadefs::num_t>(i));
    y.push_back(static_cast<datadefs::num_t>(i));
    y2.push_back(static_cast<datadefs::num_t>(49-i));
  }

  corr = math::pearsonCorrelation(x, y);
  newassert( fabs( corr - 1.0 ) < 1e-6 );

  corr = math::pearsonCorrelation(x, y2);
  newassert( fabs( corr + 1.0 ) < 1e-6 );

  // x and y are defined as const in our signature. Since we're not
  //  testing edge cases of non-trivial memory corruption, we ignore it here.

  x.clear();
  y.clear();
  y2.clear();

  corr = math::pearsonCorrelation(x, y);
  newassert( datadefs::isNAN(corr) );

  for (int i = 0; i < 50; ++i) {
    x.push_back(datadefs::NUM_NAN);
    y.push_back(datadefs::NUM_NAN);
  }

  corr = math::pearsonCorrelation(x, y);
  newassert( datadefs::isNAN(corr) );
}

void math_newtest_incrementDecrementSquaredFrequency() {

  unordered_map<cat_t,size_t> freq;
  size_t sqFreq = 0;


  // TEST: increment 1,1,1,2,2
  vector<cat_t> foo = {"1","1","1","2","2"};
  math::incrementSquaredFrequency(foo[0],freq,sqFreq);
  newassert( sqFreq == 1 );

  math::incrementSquaredFrequency(foo[1],freq,sqFreq);
  newassert( sqFreq == 4 );

  math::incrementSquaredFrequency(foo[2],freq,sqFreq);
  newassert( sqFreq == 9 );

  math::incrementSquaredFrequency(foo[3],freq,sqFreq);
  newassert( sqFreq == 10 );

  math::incrementSquaredFrequency(foo[4],freq,sqFreq);
  newassert( sqFreq == 13 );

  newassert( freq.size() == 2 );


  // TEST: decrement 2,2,1,1,1
  math::decrementSquaredFrequency(foo[4],freq,sqFreq);
  newassert( sqFreq == 10 );

  math::decrementSquaredFrequency(foo[3],freq,sqFreq);
  newassert( sqFreq == 9 );

  math::decrementSquaredFrequency(foo[2],freq,sqFreq);
  newassert( sqFreq == 4 );

  math::decrementSquaredFrequency(foo[1],freq,sqFreq);
  newassert( sqFreq == 1 );

  math::decrementSquaredFrequency(foo[0],freq,sqFreq);
  newassert( sqFreq == 0 );

  newassert( freq.size() == 0 );


  // TEST: increment 1,2,3,...,10
  for( size_t i = 1; i <= 10; ++i ) {
    math::incrementSquaredFrequency(utils::num2str(i),freq,sqFreq);
    newassert( freq.size() == i );
    newassert( sqFreq == i );
  }


  // TEST: decrement 1,2,3,...,10
  for( size_t i = 1; i <= 10; ++i ) {
    math::decrementSquaredFrequency(utils::num2str(i),freq,sqFreq);
    newassert( freq.size() == 10 - i );
    newassert( sqFreq == 10 - i );
  }


}

void math_newtest_deltaImpurity_class() {

  newassert( fabs( math::deltaImpurity_class(25900,278,9,3,25393,275) - 0.007815518677275 ) < 1e-5 );

}


#endif // MATH_NEWTEST_HPP
