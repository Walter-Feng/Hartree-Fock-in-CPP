#include "optimize.h"

#include <gsl/gsl_multimin.h>

#include "gradient/driver.h"
#include "global/error.h"
#include "util/gsl_converter.h"
#include "util/time.h"
#include "util/json.h"
#include "util/printer.h"

namespace hfincpp::optimize {

struct Optimizer {
  gradient::EnergyDriver energy_driver;
  gradient::GradientDriver gradient_driver;
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


double gradient_driver_gsl_wrapper_energy(const gsl_vector * flattened_xyz,
                                          void * param) {

  const arma::vec arma_flattened_xyz = gsl::convert_vec(flattened_xyz);
  const auto converted_optimizer = (Optimizer *) param;

  auto atoms_copy = converted_optimizer->atoms;

  const arma::mat xyz =
      arma::reshape(arma_flattened_xyz, arma::size(atoms_copy.xyz));

  atoms_copy.xyz = xyz;

  return converted_optimizer->energy_driver(atoms_copy);

}

void gradient_driver_gsl_wrapper_gradient(const gsl_vector * flattened_xyz,
                                          void * param,
                                          gsl_vector * g) {

  const arma::vec arma_flattened_xyz = gsl::convert_vec(flattened_xyz);
  const auto converted_optimizer = (Optimizer *) param;

  auto atoms_copy = converted_optimizer->atoms;

  const arma::mat xyz =
      arma::reshape(arma_flattened_xyz, arma::size(atoms_copy.xyz));


  atoms_copy.xyz = xyz;

  const auto result = converted_optimizer->gradient_driver(atoms_copy);

  auto result_pointer = gsl::convert_vec(arma::vectorise(result.second));

  gsl_vector_memcpy(g, result_pointer);
  gsl_vector_free(result_pointer);
}


void gradient_driver_gsl_wrapper(const gsl_vector * flattened_xyz,
                                 void * param,
                                 double * energy,
                                 gsl_vector * gradient) {

  const arma::vec arma_flattened_xyz = gsl::convert_vec(flattened_xyz);
  const auto converted_optimizer = (Optimizer *) param;

  auto atoms_copy = converted_optimizer->atoms;

  const arma::mat xyz =
      arma::reshape(arma_flattened_xyz, arma::size(atoms_copy.xyz));

  atoms_copy.xyz = xyz;

  const auto result = converted_optimizer->gradient_driver(atoms_copy);

  * energy = result.first;

  auto result_pointer = gsl::convert_vec(arma::vectorise(result.second));

  gsl_vector_memcpy(gradient, result_pointer);
  gsl_vector_free(result_pointer);

}


Printer<gsl_multimin_fdfminimizer *>
    generic_optimize_printer = [](
        const gsl_multimin_fdfminimizer * state,
        const double computation_time,
        const int iter,
        const int print_level,
        const bool print_header) -> int {

  int width = 18;
  int precision = 8;

  if (print_level >= 2) {
    width = 27;
    precision = 17;
  }

  int total_length = 0;

  if (print_level >= 1) {
    total_length = 6 + width * 3;

    if (print_header) {
      fmt::print("Performing optmization ...\n");
      print_separator(total_length);
      fmt::print("{:>{}}", "|Iter|", 6);
      fmt::print("{:>{}}", "Time / s |", width);
      fmt::print("{:>{}}", "Energy / a.u. |", width);
      fmt::print("{:>{}}", "Grad / a.u. |", width);
      fmt::print("\n");
      print_separator(total_length);
    }

    fmt::print("{:>{}}", iter, 6);
    const arma::rowvec combined{computation_time, state->f,
                                arma::abs(gsl::convert_vec(state->gradient)).max()};
    print(combined, width, precision);
    fmt::print("\n");
  }

  return total_length;

};

std::pair<double, geometry::Atoms>
optimize(
    const gradient::EnergyDriver & energy_driver,
    const gradient::GradientDriver & gradient_driver,
    const geometry::Atoms & atoms,
    const double initial_step_size,
    const double tolerance,
    const double gradient_tolerance,
    const size_t total_steps,
    const int print_level) {

  const auto minimizer_type = minimizer_map("bfgs");

  const size_t n_variables = atoms.xyz.n_elem;

  auto minimizer_environment =
      gsl_multimin_fdfminimizer_alloc(minimizer_type, n_variables);

  /* assigning function to minimizer object */
  gsl_multimin_function_fdf minimizer_object;
  minimizer_object.f = &gradient_driver_gsl_wrapper_energy;
  minimizer_object.df = &gradient_driver_gsl_wrapper_gradient;
  minimizer_object.fdf = &gradient_driver_gsl_wrapper;
  minimizer_object.n = n_variables;

  Optimizer param{energy_driver, gradient_driver, atoms};

  minimizer_object.params = (void *) &param;

  /* starting point */
  const arma::vec flattened = arma::vectorise(atoms.xyz);
  gsl_vector * xyz_in_gsl_vector = gsl::convert_vec(flattened);

  gsl_multimin_fdfminimizer_set(minimizer_environment,
                                &minimizer_object,
                                xyz_in_gsl_vector,
                                initial_step_size, tolerance);

  size_t iter = 0;
  int status = GSL_CONTINUE;

  const int total_length =
      generic_optimize_printer(minimizer_environment, 0.0, 0, print_level, true);

  Timer optimize_time;

  do {

    iter++;

    status = gsl_multimin_fdfminimizer_iterate(minimizer_environment);

    generic_optimize_printer(minimizer_environment,
                             optimize_time.elapsed(),
                             iter,
                             print_level,
                             false);

    if (status) {
      throw Error(gsl_strerror(status));
    }

    status = gsl_multimin_test_gradient(minimizer_environment->gradient,
                                        gradient_tolerance);

    if (status == GSL_SUCCESS) {

      if(print_level >= 1) {
        print_separator(total_length);
      }

      const arma::vec result = gsl::convert_vec(minimizer_environment->x);

      const double energy = minimizer_environment->f;

      gsl_multimin_fdfminimizer_free(minimizer_environment);
      gsl_vector_free(xyz_in_gsl_vector);

      auto atoms_copy = atoms;
      atoms_copy.xyz = arma::reshape(result, arma::size(atoms.xyz));

      return {energy, atoms_copy};
    }

  } while (status == GSL_CONTINUE && iter < total_steps);

  throw Error("fail to converge towards the solution");
}


nlohmann::json optimize(const nlohmann::json & input,
                        const geometry::Atoms & atoms,
                        const basis::Basis & basis,
                        const std::string & method) {
  const auto gradient_driver = gradient::driver(input, atoms, basis, method);
  const auto energy_driver = gradient::energy_driver(input, atoms, basis, method);

  nlohmann::json json_output;

  const int print_level = input["print_level"];

  const auto optimized_geometry = optimize(energy_driver, gradient_driver, atoms);

  util::put(json_output, "geometry", optimized_geometry.second.xyz);
  util::put(json_output, "energy", optimized_geometry.first);

  if(print_level >= 1) {
    int width = 18;
    int precision = 8;

    fmt::print("Optimized geometry:\n");
    if (print_level >= 2) {
      width = 27;
      precision = 17;
    }

    int total_length = 6 + width * 3;

    print_separator(total_length);
    fmt::print("{:>{}}", "|Atom|", 6);
    fmt::print("{:>{}}", "X / Angstrom |", width);
    fmt::print("{:>{}}", "Y / Angstrom |", width);
    fmt::print("{:>{}}", "Z / Angstrom |", width);
    fmt::print("\n");
    print_separator(total_length);

    for(int atom=0; atom < atoms.n_atoms(); atom++) {
      fmt::print("{:>{}}", atoms.symbols[atom], 6);
      print(arma::rowvec(optimized_geometry.second.xyz.col(atom).t()), width, precision);
      fmt::print("\n");
    }
    print_separator(total_length);
  }


  return json_output;

}


}