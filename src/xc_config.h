// XC configuration for unified SCF. HF is modeled via exchange="HF", correlation="NONE", ax=1.0
#pragma once

#include <string>

struct XCConfig {
  // Exchange functional identifier: "HF" for exact exchange or libxc id/name
  std::string exchange = "HF";
  // Correlation functional identifier: "NONE" or libxc id/name
  std::string correlation = "NONE";
  // Fraction of exact exchange (hybrids). 1.0 for HF, 0.0 for pure LDA/GGA.
  double exact_exchange_fraction = 1.0;
};
