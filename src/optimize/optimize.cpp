#include "optimize.h"

#include <gsl/gsl_multimin.h>

#include "util/gsl_converter.h"

namespace hfincpp::optimize {

struct Optimizer {
  gradient::GradientDriver driver;
  geometry::Atoms atoms;
};

gsl_multimin_fdfminimizer_type * minimizer_map(const std::string type) {
  if(type == "steepest_descent") {
    return const_cast<gsl_multimin_fdfminimizer_type *>(gsl_multimin_fdfminimizer_steepest_descent);
  } else if(type == "conjugate_pr") {
    return const_cast<gsl_multimin_fdfminimizer_type *>(gsl_multimin_fdfminimizer_conjugate_pr);
  } else if(type == "conjugate_fr") {
    return const_cast<gsl_multimin_fdfminimizer_type *>(gsl_multimin_fdfminimizer_conjugate_fr);
  } else if(type == "bfgs") {
    return const_cast<gsl_multimin_fdfminimizer_type *>(gsl_multimin_fdfminimizer_vector_bfgs);
  } else if(type == "bfgs2") {
    return const_cast<gsl_multimin_fdfminimizer_type *>(gsl_multimin_fdfminimizer_vector_bfgs2);
  } else {
    throw Error("minimizer " + type + " is not implemented");
  }

}

void gradient_driver_gsl_wrapper(const gsl_vector * flattened_xyz,
                                 void * param,
                                 double * energy,
                                 gsl_vector * gradient) {

  const arma::vec arma_flattened_xyz = gsl::convert_vec(flattened_xyz);
  const auto converted_optimizer = (Optimizer *) param;

  auto atoms_copy = converted_param->atoms;

  const arma::mat xyz =
      arma::reshape(arma_flattened_xyz, arma::size(atoms_copy.xyz));

  atoms_copy.xyz = xyz;

  const auto result = converted_optimizer->driver(atoms_copy);

  * energy = result.first;

  const auto result_pointer = gsl::convert_vec(result);

  gsl_vector_memcpy(g, result_pointer);

  gsl_vector_free(result_pointer);

}

geometry::Atoms optimize(const gradient::GradientDriver & gradient_driver,
                         const geometry::Atoms & atoms,
                         const std::string optimization_method) {

}

}