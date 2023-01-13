#include "gradient.h"

#include "hf/rhf.h"
#include "util/printer.h"

namespace hfincpp::gradient {

nlohmann::json driver(const nlohmann::json & input,
                      const geometry::Atoms & atoms,
                      const basis::Basis & basis,
                      const std::string method) {
  if(method == "rhf") {
    const auto gradient_driver = hf::gradient::driver(input, atoms, basis);

    const auto result = gradient_driver(atoms);

    nlohmann::json json_output;
    util::put(json_output, "energy", result.first);
    util::put(json_output, "gradient", result.second);

    const int print_level = input["print_level"];
    if(print_level >= 1) {
      int width = 18;
      int precision = 8;

      if (print_level >= 2) {
        width = 27;
        precision = 17;
      }

      int total_length = 6 + width * 3;

      for (int i = 0; i < total_length; i++) {
        fmt::print("=");
      }
      fmt::print("\n");

      fmt::print("{:>{}}", "|Atom|", 6);
      fmt::print("{:>{}}", "X / (Ha / bohr) |", width);
      fmt::print("{:>{}}", "Y / (Ha / bohr) |", width);
      fmt::print("{:>{}}", "Z / (Ha / bohr) |", width);
      fmt::print("\n");
      for (int i = 0; i < total_length; i++) {
        fmt::print("=");
      }
      fmt::print("\n");

      for(int atom=0; atom < atoms.n_atoms(); atom++) {
        fmt::print("{:>{}}", atoms.symbols[atom], 6);
        print(arma::rowvec(result.second.col(atom).t()), width, precision);
        fmt::print("\n");
      }
      for (int i = 0; i < total_length; i++) {
        fmt::print("=");
      }
      fmt::print("\n");
    }


    return json_output;
  } else {
    throw Error("method not recognized");
  }
}

}