#ifndef EXCEPTIONS_HPP
#define EXCEPTIONS_HPP

#include <exception>
#include <string>
#include <sstream>

#include "errno.hpp"

class RFACE_EXCEPTION : public exception {

public:

  RFACE_EXCEPTION(int errno, const string note = ""):errno_(errno) {
    stringstream ss;
    ss << "ERRNO (" << errno_ << "): ";
    switch ( errno_ ) {
    case ERRNO_INVALID_ARGUMENT:
      ss << "invalid command-line argument.";
      break;
    case ERRNO_INVALID_VALUE:
      ss << "invalid command-line value.";
      break;
    case ERRNO_ILLEGAL_MEMORY_ACCESS:
      ss << "illegal memory access.";
      break;
    default:
      ss << "unknown exception!";
      break;
    }
    ss << " " << note;
    msg_ = ss.str();
  }

  ~RFACE_EXCEPTION() throw() {}

  virtual const char* what() const throw() {
    return msg_.c_str();
  }
  
private:
  std::string msg_;
  int errno_;
  
};

#endif
