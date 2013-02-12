#ifndef DATADEFS_NEWTEST_HPP
#define DATADEFS_NEWTEST_HPP

#include <cstdlib>
#include <limits>

#include "newtest.hpp"
#include "datadefs.hpp"

void datadefs_newtest_isnan();

void datadefs_newtest() {

  newtest( "Testing interpretation of missing values", &datadefs_newtest_isnan );

}

void datadefs_newtest_isnan() {
 
  newassert(datadefs::isNAN_STR("na"));
  newassert(datadefs::isNAN_STR("NA"));
  newassert(datadefs::isNAN_STR("nan"));
  newassert(datadefs::isNAN_STR("NaN"));
  newassert(datadefs::isNAN_STR("NAN"));
  newassert(datadefs::isNAN_STR("?"));
  
  newassert(!datadefs::isNAN_STR("2"));
  newassert(!datadefs::isNAN_STR("@data"));
  newassert(!datadefs::isNAN_STR("NAte"));
  newassert(!datadefs::isNAN_STR("name"));
  newassert(!datadefs::isNAN_STR("na?"));

  newassert(numeric_limits<datadefs::num_t>::has_quiet_NaN);
  newassert(datadefs::isNAN(datadefs::NUM_NAN));
  
  if (numeric_limits<datadefs::num_t>::has_infinity) {
    newassert(!datadefs::isNAN(numeric_limits<datadefs::num_t>::infinity()));
    newassert(!datadefs::isNAN(-numeric_limits<datadefs::num_t>::infinity()));
  }
  
  newassert(!datadefs::isNAN((datadefs::num_t)0.0));
  newassert(!datadefs::isNAN((datadefs::num_t)-1.0));
  newassert(!datadefs::isNAN((datadefs::num_t)1.0));
  
}

#endif
