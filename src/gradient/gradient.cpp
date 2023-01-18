#include "gradient.h"

#include "hf/rhf.h"
#include "util/printer.h"

namespace hfincpp::gradient {

nlohmann::json gradient(const nlohmann::json & input,
                        const geometry::Atoms & atoms,
                        const basis::Basis & basis,
                        const std::string method) {

  const int print_level = input["print_level"];

  const auto gradient_driver = driver(input, atoms, basis, method);

  if(print_level >= 1) {
    fmt::print("Calculating gradient ...\n");
  }

  const auto result = gradient_driver(atoms);

  nlohmann::json json_output;
  util::put(json_output, "energy", result.first);
  util::put(json_output, "gradient", result.second);


  if(print_level >= 1) {
    int width = 18;
    int precision = 8;

    if (print_level >= 2) {
      width = 27;
      precision = 17;
    }

    int total_length = 6 + width * 3;

    if(print_level >= 1) {
      print_separator(total_length);
    }

    fmt::print("{:>{}}", "|Atom|", 6);
    fmt::print("{:>{}}", "X / (Ha / bohr) |", width);
    fmt::print("{:>{}}", "Y / (Ha / bohr) |", width);
    fmt::print("{:>{}}", "Z / (Ha / bohr) |", width);
    fmt::print("\n");
    print_separator(total_length);
    for(int atom=0; atom < atoms.n_atoms(); atom++) {
      fmt::print("{:>{}}", atoms.symbols[atom], 6);
      print(arma::rowvec(result.second.col(atom).t()), width, precision);
      fmt::print("\n");
    }
    if(print_level >= 1) {
      print_separator(total_length);
    }
  }


  return json_output;

}

}