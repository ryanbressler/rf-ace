#ifndef ARGPARSE_HPP
#define ARGPARSE_HPP

#include <cassert>
#include <iostream>
#include <map>
#include <sstream>
#include <string.h>
#include <vector>
#include "errno.hpp"

using namespace std;

/**
 * Generic argument parser that allows options to be checked
 * speculatively. Currently uses a map representation of the argument tree. 
 */
class ArgParse {

public:
  ArgParse(const int argc, char* const argv[]) {
    if (argc < 1) {
      throw ERRNO_INVALID_ARGUMENT;
    }

    string currArg = "";
    for (int i = 0; i < argc; ++i) {
      
      try {
        // !! Correctness Note: this strategy may attempt to dereference
        // !!  corrupt memory. Thus, these runtime checks, while better than an
        // !!  outright crash, may imply security vulnerabilities in dependent
        // !!  code. Beware!
        if (argv[i] == NULL || argv[i][0] == 0x0) {
          throw ERRNO_INVALID_ARGUMENT;
        } else {
          if (!currArg.empty()) {
            mappedArgs[currArg] = string(argv[i]);
            currArg = string("");
            //continue;
          } 
          switch(argv[i][0]) {
          case '-':
          case '+':
            if (argv[i][1] == '\0') { break; }
            if (argv[i][1] == '-' ||
                argv[i][1] == '+') {
              
              stringstream argSS;
              bool containsEquals = false;
              size_t idx;
              size_t len = strlen(argv[i]);
              for (idx = 2; idx < len; ++idx) {
                if (argv[i][idx] == '=') {
                  containsEquals = true;
                  break;
                }
                argSS << (char)argv[i][idx];
              }
              
              currArg = argSS.str();
              if (containsEquals) {
                stringstream valSS;
                for (++idx; idx < len; ++idx) {
                  valSS << (char)argv[i][idx];
                }
                mappedArgs[currArg] = valSS.str();
                currArg = string("");
              }
              
              
            } else {
              size_t len = strlen(argv[i]);
              for (size_t idx = 1; idx < len; ++idx) {
                char arg[] = {argv[i][idx], '\0'};
                mappedArgs[string(arg)] = string("");
                currArg = string(arg);
              }
            }
            
            break;
          default:
            extraArgs.push_back(string(argv[i]));
            break;
          }
        }
      } catch (...) {
        //assert(false); // Check if argv was corrupt or overstepped. Implies a
                         //  major FIXME if this is hit. (Disabled in lieu of
                         //  runtime checks during testing)
        
        throw ERRNO_ILLEGAL_MEMORY_ACCESS;  // Perform a safer runtime check
                                            //  that should never be hit by
                                            //  correct code.
      }
    }
    if (!currArg.empty()) {
      mappedArgs[currArg] = string("");
      currArg = string("");
    }

    /*
    cout << "Mapped args:" << endl;
    for (map<string,string>::iterator it = mappedArgs.begin(); it != mappedArgs.end(); ++it) {
      cout << (*it).first << "->" << (*it).second << endl;
    }

    cout << "Extra args:" << endl;
    for (int i = 0; i < extraArgs.size(); ++i) {
      cout << extraArgs[i] << endl;
      }*/
  }
  ~ArgParse() {}

  /**
   * Queries the backend map for the current argument-value pair. Extra
   *  arguments passed positionally are not yet supported.
   *
   *  Contractual guarantees:
   *
   *  + The behavior of the extraction operator (>>) will be used with input
   *     type T. This may cause unexpected results if your type explicitly
   *     specifies an append instead of an overwrite for this
   *     operation. Declare the contents of returnVal carefully or redefine
   *     your type for these cases.
   *
   *  + Certain types specifiable for T may cause memory access violations that
   *     are difficult to debug. For example, specifying char* may throw a
   *     memory access violation at 'SOptarg >> returnVal'. It is expected and
   *     indeed required that your types be well-defined before passing them
   *     to this method.
   *
   *  + Attempting to pass a NULL pointer for any input value will not
   *     work. Don't do it.
   *
   *  + Duplicate arguments will prefer the last long specification over the
   *     last short specification of that same argument.
   *
   *  + Arguments are case-sensitive.
   *
   *  Sets returnVal and returns true if an argument was found; false
   *  otherwise. !! TODO return a unified status code instead
   */
  template <typename T> bool getArgument(const char* shortName, const char* longName, T& returnVal) {

    assert(shortName != NULL);
    assert(longName != NULL);
    assert(strlen(shortName) == 1);
    assert(*longName != 0);
    
    map<string,string>::iterator it = mappedArgs.find(longName);
    if (it == mappedArgs.end()) {
        it = mappedArgs.find(shortName);
    }

    if (it != mappedArgs.end()) {
      string found = (*it).second;
      if (found.empty()) {
        throw ERRNO_INVALID_VALUE;
      }
      stringstream ss(found);
      ss >> returnVal;

      if (ss.fail() || !ss.eof()) {
        throw ERRNO_INVALID_VALUE;
      }
      return true;
    }
    
    return false;
  }

  template <typename T> bool getArgument(string& shortName, string& longName, T& returnVal) {
    return getArgument<T>(shortName.c_str(), longName.c_str(), returnVal);
  }

  /**
   * Queries the backend map for the presence of a flag. This can also be used
   *  to check for the presence of an argument-value pair or non-presence of a
   *  value, but abusing this functionality is not recommended.
   *
   *  Contractual guarantees:
   *
   *  + Attempting to pass a NULL pointer for any input value will not
   *     work. Don't do it.
   *
   *  + Duplicate flags are assumed to be one instance of the set
   *     flag. Conflicting, non-duplicate flags are your problem. 
   *
   *  + Flags are case-sensitive.
   *
   */
  bool getFlag(const char* shortName, const char* longName, bool& returnVal) {
    
    assert(shortName != NULL);
    assert(longName != NULL);
    assert(strlen(shortName) == 1);
    assert(*longName != 0);
    
    map<string,string>::iterator it = mappedArgs.find(longName);
    if (it == mappedArgs.end()) {
        it = mappedArgs.find(shortName);
    }

    if (it != mappedArgs.end()) {
      returnVal = true;
      return true;
    }
    
    return false;
  }

  bool getFlag(string& shortName, string& longName, bool& returnVal) {
    return getFlag(shortName.c_str(), longName.c_str(), returnVal);
  }

  
private:
  map<string,string> mappedArgs;
  vector<string> extraArgs;
};

#endif
