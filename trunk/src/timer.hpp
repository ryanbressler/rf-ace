#ifndef TIMER_HPP
#define TIMER_HPP

#include <cstdlib>
#include <vector>
#include <string>
#include <cassert>
#include <ctime>
#include <cmath>
#include <map>

#include "datadefs.hpp"

using namespace std;

class Timer {

public:

  Timer() {}
  ~Timer() {}

  void tic(const string& objName) { 
    name2idx_.insert( pair<string,size_t>(objName,timedObjects_.size()) );
    timedObjects_.push_back( TimedObject(objName) );
  }

  void toc(const string& objName) {

    map<string,size_t>::const_iterator it( name2idx_.find(objName) );

    if ( it == name2idx_.end() ) {
      cerr << "Cannot stop timing '" << objName << "', since it was never started!" << endl;
      exit(1);
    }

    size_t idx = name2idx_[objName];
    timedObjects_[idx].timeDiff = time(0) - timedObjects_[idx].startTime;
    timedObjects_[idx].clockDiff = clock() - timedObjects_[idx].startClocks;
    if ( timedObjects_[idx].timeDiff > 0 ) {
      timedObjects_[idx].boost = static_cast<clock_t>(round(1.0 * timedObjects_[idx].clockDiff / ( CLOCKS_PER_SEC * timedObjects_[idx].timeDiff )));
    }
    timedObjects_[idx].isRunning = false;
  }

  void print() {
    cout << "Execution time breakdown:" << endl;
    for ( size_t i = 0; i < timedObjects_.size(); ++i ) {
      timedObjects_[i].print();
    }
    cout << endl;
  }

private:

  struct TimedObject {
    string  name;
    clock_t startTime;
    clock_t timeDiff;
    clock_t startClocks;
    clock_t clockDiff;
    clock_t boost;
    bool    isRunning;
    TimedObject(const string& newName): name(newName),startTime(time(0)),startClocks(clock()),boost(1),isRunning(true) {}
    void print() {
      if ( !isRunning ) {
	cout << name << "  " << timeDiff << " seconds (" << boost << "x)" << endl;
      } else {
	cout << name << " is still running!" << endl;
      }
    }
  };
  
  map<string,size_t> name2idx_;

  vector<TimedObject> timedObjects_;

};


#endif
