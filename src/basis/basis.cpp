#include "basis.h"

#include <boost/math/special_functions/factorials.hpp>
#include "geometry/periodic_table.h"
#include "util/json.h"
#include "global/error.h"

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

  this->basis_name = basis_name;

  if (f.is_open()) {

    const nlohmann::json data = nlohmann::json::parse(f);

    Basis basis;

    std::vector<arma::uword> atom_indices_in_std_vector;
    std::vector<arma::uword> shell_indices_in_std_vector;
    const int n_atoms = atoms.n_atoms();
    for (int i_atom = 0; i_atom < n_atoms; i_atom++) {

      const int atomic_numbers = atoms.atomic_numbers[i_atom];
      const arma::vec3 center = atoms.xyz.col(i_atom);
      const nlohmann::json atom_json =
          data["elements"][std::to_string(atomic_numbers)];

      for (auto ith_layer: atom_json["electron_shells"]) {

        const auto ith_layer_coefficients = ith_layer["coefficients"];

        const auto exponents_in_strings = ith_layer["exponents"];

        for (size_t ith_shell = 0;
             ith_shell < ith_layer_coefficients.size(); ith_shell++) {

          const auto angular_momentum = ith_layer["angular_momentum"];

          /* There are two formats of basis sets, one being Pople
           * and not Pople (for example, cc-pvdz).
           * The formats are just different, with Pople ones stacking SPDF..
           * into one set of exponent / coefficient pairs for the same
           * principal quantum number (n), while others stack orbitals with
           * different n for a specific s/p/d/f orbital.
           *
           * Go check data/basis_set_exchange if I failed to explain
           * it explicitly.
           * */
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

          const size_t n_functions = coefficients_in_strings.size();

          assert(exponents_in_strings.size() == n_functions);

          arma::vec coefficients(n_functions);
          arma::vec exponents(n_functions);

          for (size_t i = 0; i < n_functions; i++) {
            const std::string coefficient = coefficients_in_strings[i];
            const std::string exponent = exponents_in_strings[i];
            coefficients(i) = std::stod(coefficient);
            exponents(i) = std::stod(exponent);
          }

          const arma::uvec non_zero_coefficients = arma::find(coefficients);

          const GTOShell shell{center,
                               i_angular_momentum,
                               exponents(non_zero_coefficients)};

          shells.push_back(shell);

          std::string shell_label = std::to_string(i_atom + 1) + atoms.symbols[i_atom];
          shell_label += angular_string;
          shell_labels.push_back(shell_label);

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
          for (int ax = i_angular_momentum; ax >= 0; ax--) {

            const std::string ax_string = ax ? "x^" + std::to_string(ax) : "";

            for (int ay = i_angular_momentum - ax; ay >= 0; ay--) {

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

              std::string label = std::to_string(i_atom + 1) + atoms.symbols[i_atom] + " ";
              label += angular_string;
              label += ax_string;
              label += ay_string;
              label += az_string;

              const GTOFunction function{center, angular,
                                         exponents(non_zero_coefficients),
                                         coefficients(non_zero_coefficients)
                                         % normalization_constant(non_zero_coefficients)};

              functions.push_back(function);
              function_labels.push_back(label);
              shell_indices_in_std_vector.push_back(shells.size()-1);
              atom_indices_in_std_vector.push_back(i_atom);
            }
          }
        }
      }
    }

    atomic_numbers = atoms.atomic_numbers;
    atom_symbols = atoms.symbols;
    atom_indices = arma::uvec{atom_indices_in_std_vector};
    shell_indices = arma::uvec{shell_indices_in_std_vector};
    to_conventional_function_indexing =
        arma::linspace<arma::uvec>(0, this->n_functions() - 1, 1);

    *this = this->sort_by_angular_momentum();

  } else {
    throw Error("The basis is not found");
  }
}

int Basis::n_atoms() const {
  return atomic_numbers.n_elem;
}

int Basis::n_shells() const {
  return shells.size();
}

int Basis::n_functions() const {
  return functions.size();
}

arma::uvec Basis::on_atom(const arma::uword atom_index) const {
  return arma::find(atom_indices == atom_index);
}

std::vector<arma::uvec> Basis::on_atoms() const {
  std::vector<arma::uvec> result(n_atoms());
for(int i=0; i<n_atoms(); i++) {
    result[i] = on_atom(i);
  }
  return result;
}

Basis Basis::sort_by_angular_momentum() const {
  const auto n_functions = functions.size();
  const auto n_shells = shells.size();
  arma::uvec angular_momentum(n_functions);
  arma::uvec angular_momentum_in_shells(n_shells);

  for(size_t i=0; i<n_functions; i++) {
    angular_momentum(i) = arma::sum(functions[i].angular);
  }

  for(size_t i=0; i<n_shells; i++) {
    angular_momentum_in_shells(i) = arma::sum(shells[i].angular);
  }

  Basis new_basis;
  arma::uvec function_indexing;

  for(arma::uword l=0; l<=arma::max(angular_momentum); l++) {
    const arma::uvec find = arma::find(angular_momentum == l);
    function_indexing = arma::join_vert(function_indexing, find);
    const arma::uvec find_shell = arma::find(angular_momentum_in_shells == l);

    assert(find.is_sorted());
    assert(find_shell.is_sorted());

    for(arma::uword i=0; i<find.n_elem; i++) {
      new_basis.functions.push_back(functions[find(i)]);
      new_basis.function_labels.push_back(function_labels[find(i)]);
    }

    for(arma::uword i=0; i<find_shell.n_elem; i++) {
      new_basis.shells.push_back(shells[find_shell(i)]);
      new_basis.shell_labels.push_back(shell_labels[find_shell(i)]);
    }
  }

  assert(function_indexing.n_elem == atom_indices.n_elem);
  new_basis.atomic_numbers = atomic_numbers;
  new_basis.atom_symbols = atom_symbols;
  new_basis.atom_indices = atom_indices(function_indexing);
  new_basis.shell_indices = shell_indices(function_indexing);
  new_basis.basis_name = basis_name;
  new_basis.to_conventional_function_indexing =
      to_conventional_function_indexing(function_indexing);
  return new_basis;
}
}
}