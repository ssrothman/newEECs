#ifndef SRothman_EECs_Warn_h
#define SRothman_EECs_Warn_h

#ifdef CMSSW_GIT_HASH
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#else
#include <iostream>
#endif

namespace eec::detail {

  // CMSSW_GIT_HASH is defined by scram builds only.
  inline void warn(const char* msg) {
#ifdef CMSSW_GIT_HASH
    edm::LogWarning("EECs") << msg;
#else
    std::cerr << "eec warning: " << msg << '\n';
#endif
  }

}  // namespace eec::detail

#endif
