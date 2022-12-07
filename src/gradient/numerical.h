#ifndef GRADIENT_NUMERICAL_H
#define GRADIENT_NUMERICAL_H

#include <armadillo>
#include <functional>

#include "geometry/geometry.h"

namespace hfincpp::gradient {

template<typename T>
std::vector<T>
numerical(const std::function<T(const arma::vec &)> & functor,
          const arma::vec & reference_point,
          const double step_size = 1e-4) {
  std::vector<T> result(reference_point.n_elem);

  for(size_t i=0; i<result.size(); i++) {
    arma::vec forward = reference_point;
    arma::vec backward = reference_point;
    forward(i) += step_size;
    backward(i) -= step_size;

    result[i] = (functor(forward) - functor(backward)) / 2.0 / step_size;
  }

  return result;
}

template<typename T>
std::vector<T>
numerical(const std::function<T(const geometry::Atoms &)> & functor,
          const geometry::Atoms & reference,
          const double step_size = 1e-4) {
  const arma::mat ref_xyz = reference.xyz;

  const std::function<T(const arma::vec &)> wrapped_functor =
      [functor, reference](const arma::vec & flattened_xyz) -> T {
    geometry::Atoms atoms_copy = reference;
    const arma::mat reshaped_xyz = arma::reshape(flattened_xyz, 3, atoms_copy.n_atoms());
    atoms_copy.xyz = reshaped_xyz;
    return functor(atoms_copy);
  };

  return numerical(wrapped_functor, arma::vectorise(ref_xyz), step_size);
}

}

#endif //GRADIENT_NUMERICAL_H
