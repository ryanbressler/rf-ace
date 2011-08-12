#ifndef ARGPARSE_HPP
#define ARGPARSE_HPP

#include <getopt.h>
#include <cassert>
#include <iostream>
#include <sstream>
#include <string.h>
#include "errno.hpp"

using namespace std;

/**
 * Generic argument parser that allows options to be checked
 *  speculatively. Currently wraps GNU getopt_long, which is completely
 *  thread-unsafe. Let the maintainer beware.
 */
class ArgParse {

public:
  ArgParse(const int argc, char* const argv[]) {
    if (argc < 1) {
      throw ERRNO_INVALID_ARGUMENT;
    }
    
    for (int i = 0; i < argc; ++i) {
      
      try {
        // !! Correctness Note: this strategy may attempt to dereference
        // !!  corrupt memory. Thus, these runtime checks, while better than an
        // !!  outright crash, may imply security vulnerabilities in dependent
        // !!  code. Beware!
        if (argv[i] == NULL) {
          throw ERRNO_INVALID_ARGUMENT;
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
    argc_ = argc;
    argv_ = argv;
    opterr = 0; // Suppress invalid print messages from getopt, which trigger
                //  spuriously from this class when long options are issued.
  }
  ~ArgParse() {}

  /**
   * Shims the local version of ArgParse to GNU getopt_long. Changes the input
   *  assumptions for long arguments slightly, but makes the whole system more
   *  robust as a result.
   *
   *  Note that due to interface imparity, each pass of getArgument takes linear
   *  time with respect to the original command line. It also incurs the
   *  overhead of a stringstream should an option be found.
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
   *  Sets returnVal and returns true if an argument was found; false
   *  otherwise. !! TODO return a unified status code instead
   */
  template <typename T> bool getArgument(const char* shortName, const char* longName, T& returnVal) {

    assert(shortName != NULL);
    assert(longName != NULL);
    assert(strlen(shortName) == 1);
    assert(*longName != 0);
    
    const char cShortName = *shortName;
    const struct option long_options[] = {
      {longName, 1, NULL, cShortName}
    };
    const char opts[] = {cShortName, ':', '\0'};

    // Copy the contents of argv into a temporary, rearrangeable variable
    char** targv = new char* [argc_];
    for (int i = 0; i < argc_; ++i) {
      targv[i] = argv_[i];
    }
    
    int option_index = 0;
    int c = 0;
    bool rst = false;
    while (c != -1) {
      c = getopt_long (argc_, targv, opts, long_options, &option_index);
      if (c == cShortName) {
        stringstream sOptarg(optarg); // Streamify the arg received from getopt
        sOptarg >> returnVal;
        rst = true;
        break;
      } 
    }
    optind = 0; // Reset the getopt parsing index to 0
    delete[] targv; // Free our temporary argv
    return rst;
  }

  template <typename T> bool getArgument(string& shortName, string& longName, T& returnVal) {
    return getArgument<T>(shortName.c_str(), longName.c_str(), returnVal);
  }

  /**
   * Shims the local version of ArgParse to GNU getopt_long. Changes the input
   *  assumptions for long arguments slightly, but makes the whole system more
   *  robust as a result.
   *
   *  Note that due to interface imparity, each pass of getArgument takes linear
   *  time with respect to the original command line.
   *
   *  Tests for the presence of a flag. Unlike the above, returnVal must be of
   *  type bool. !! TODO return a unified status code instead
   *
   *  Contractual guarantees:
   *
   *  + Attempting to pass a NULL pointer for any input value will not
   *     work. Don't do it.
   *
   */
  bool getFlag(const char* shortName, const char* longName, bool& returnVal) {
    
    assert(strlen(shortName) == 1);
    assert(*longName != 0);

    const char cShortName = *shortName;
    const struct option long_options[] = {
      {longName, 0, NULL, cShortName}
    };

    // Copy the contents of argv into a temporary, rearrangeable variable
    char** targv = new char* [argc_];
    for (int i = 0; i < argc_; ++i) {
      targv[i] = argv_[i];
    }

    int option_index = 0;
    int c = 0;
    bool rst = false;
    returnVal = false;
    while (c != -1) {
      c = getopt_long (argc_, targv, shortName, long_options, &option_index);
      if (c == cShortName) {
        returnVal = true;
        rst = true;
        break;
      }
    }

    // Set returnVal to false to signal the non-presence of this flag
    if (!rst) {
      returnVal = false;
    }

    optind = 0; // Reset the getopt parsing index to 0
    delete[] targv; // Free our temporary argv
    return rst;
  }

  bool getFlag(string& shortName, string& longName, bool& returnVal) {
    return getFlag(shortName.c_str(), longName.c_str(), returnVal);
  }

  
private:
  int argc_;
  char* const* argv_;
};

#endif
