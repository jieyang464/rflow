#pragma once

#include "function.h"
#include "matrix.h"

namespace commutator {

enum class ConjugationMethod {
  Taylor,
  BCH,
  Injected
};

struct Options;
using ConjugationFn = scfcpp::Function<void(const scfcpp::DenseMatrix&,
                                            const scfcpp::DenseMatrix&,
                                            const Options&,
                                            scfcpp::DenseMatrix&)>;

struct Options {
  double alpha = 0.1;
  int exp_order = 4;
  int bch_order = 4;
  bool do_purify = false;
  int purify_steps = 1;
  ConjugationMethod method = ConjugationMethod::Taylor;
  ConjugationFn conjugator;
};

// Perform one RHF AO S-metric commutator step on total density P (where spin density is P/2)
// P is updated in-place.
void update_step_rhf(const scfcpp::DenseMatrix& F,
                     const scfcpp::DenseMatrix& S,
                     const Options& opts,
                     scfcpp::DenseMatrix& P);

// Perform one UHF AO S-metric commutator step on spin densities Pa and Pb separately
// Pa and Pb are updated in-place.
void update_step_uhf(const scfcpp::DenseMatrix& Fa,
                     const scfcpp::DenseMatrix& Fb,
                     const scfcpp::DenseMatrix& S,
                     const Options& opts,
                     scfcpp::DenseMatrix& Pa,
                     scfcpp::DenseMatrix& Pb);

}
