#ifndef SCF_SETUP_H
#define SCF_SETUP_H

#include "scf.h"

#include "geometry/geometry.h"
#include "basis/basis.h"

namespace hfincpp::scf {

template<typename T>
struct Setup {
  geometry::Atoms atoms;
  basis::Basis basis;
  scf::OverlapMatrix<T> overlap;

  int print_level;
};

template<typename T>
struct Result {
  Setup<T> setup;
  scf::Result<T> result;
};

}

#endif //SCF_SETUP_H
