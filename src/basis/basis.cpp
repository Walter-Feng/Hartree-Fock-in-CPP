#include "basis.h"

#include <boost/math/special_functions/factorials.hpp>
#include "geometry/periodic_table.h"
#include "util/json.h"
#include "util/error.h"

namespace hfincpp {
namespace basis {

Basis::Basis() {
  functions = {};
  atom_symbols = {};
  atomic_numbers = arma::uvec{};
  function_labels = {};
}

Basis::Basis(const geometry::Atoms & atoms,
             const std::string & basis_name) {

  const std::string data_path = BASIS_PATH;

  std::ifstream f(data_path + basis_name + ".0.json");

  if (f.is_open()) {

    const nlohmann::json data = nlohmann::json::parse(f);

    Basis basis;

    const int n_atoms = atoms.n_atoms();
    for (int i_atom = 0; i_atom < n_atoms; i_atom++) {

      const int atomic_numbers = atoms.atomic_numbers[i_atom];
      const arma::vec3 center = atoms.xyz.col(i_atom);
      const nlohmann::json atom_json =
          data["elements"][std::to_string(atomic_numbers)];

      for (auto ith_layer: atom_json["electron_shells"]) {

        const auto ith_layer_coefficients = ith_layer["coefficients"];

        const auto exponents_in_strings = ith_layer["exponents"];

        for (int ith_shell = 0;
             ith_shell < ith_layer_coefficients.size(); ith_shell++) {

          const auto angular_momentum = ith_layer["angular_momentum"];
          const bool is_pople_basis_set = angular_momentum.size() > 1;

          const int i_angular_momentum =
              is_pople_basis_set ?
              angular_momentum[ith_shell].get<int>() :
              angular_momentum[0].get<int>();

          std::string angular_string;
          switch (i_angular_momentum) {
            case 0:
              angular_string = "s";
              break;
            case 1:
              angular_string = "p";
              break;
            case 2:
              angular_string = "d";
              break;
            case 3:
              angular_string = "f";
              break;
            case 4:
              angular_string = "g";
              break;
            default:
              throw Error("unknown angular momentum number encountered");
          }

          const auto coefficients_in_strings =
              ith_layer_coefficients[ith_shell];

          const int n_functions = coefficients_in_strings.size();

          assert(exponents_in_strings.size() == n_functions);

          arma::vec coefficients(n_functions);
          arma::vec exponents(n_functions);

          for (int i = 0; i < n_functions; i++) {
            const std::string coefficient = coefficients_in_strings[i];
            const std::string exponent = exponents_in_strings[i];
            coefficients(i) = std::stod(coefficient);
            exponents(i) = std::stod(exponent);
          }
          /*
           * Hey, think of this as follows:
           * You would like to split some balls, total number being
           * i_angular_momentum, by three boxes, so that it becomes
           *       o   o   o      <--- balls
           *  x |  y |   z        <--- boxes
           *  then the possible positions of seperators, which are both identical,
           *  are (i_angular_momentum + 1) * (i_angular_momentum + 2) / 2.
           *         ^ insert first separator   ^ second seperator     ^ identical
           *  By the way, it was yz^2 orbital among f orbitals.
           *  You can also do, you first take arbitrary balls for x,
           *  then you decide how to split between y and z.
           *  This will give you 1 + 2 + ... + (i_angular_momentum + 1).
           *  Yeah these two match.
           */
          for (int ax = 0; ax <= i_angular_momentum; ax++) {

            const std::string ax_string = ax ? "x^" + std::to_string(ax) : "";

            for (int ay = 0; ay <= i_angular_momentum - ax; ay++) {

              const std::string ay_string = ay ? "y^" + std::to_string(ay) : "";

              const int az = i_angular_momentum - ax - ay;
              const std::string az_string = az ? "z^" + std::to_string(az) : "";

              const arma::Col<int>::fixed<3> angular = {ax, ay, az};

              const double scalar_factor = std::sqrt(
                  boost::math::factorial<double>(ax)
                  * boost::math::factorial<double>(ay)
                  * boost::math::factorial<double>(az)
                  / (boost::math::factorial<double>(ax * 2)
                     * boost::math::factorial<double>(ay * 2)
                     * boost::math::factorial<double>(az * 2))
              );

              const arma::vec normalization_constant =
                  arma::pow(2.0 * exponents / M_PI, 0.75)
                  % arma::pow(8.0 * exponents, i_angular_momentum / 2.0) *
                  scalar_factor;

              std::string label = atoms.symbols[i_atom] + " ";
              label += angular_string;
              label += ax_string;
              label += ay_string;
              label += az_string;

              const GTOFunction function{center, angular, exponents,
                                         coefficients % normalization_constant};

              functions.push_back(function);
              function_labels.push_back(label);
            }
          }


        }
      }

    }

    atomic_numbers = atoms.atomic_numbers;
    atom_symbols = atoms.symbols;

  } else {
    throw Error("The basis is not found");
  }
}

Basis::Basis(const Basis & basis) {
  functions = basis.functions;
  atom_symbols = basis.atom_symbols;
  atomic_numbers = basis.atomic_numbers;
}

int Basis::n_atoms() const {
  return atomic_numbers.n_elem;
}

int Basis::n_functions() const {
  return functions.size();
}

}
}