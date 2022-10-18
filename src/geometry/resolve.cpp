#include "resolve.h"

namespace hfincpp {
namespace util {

template<>
geometry::Atoms resolve(const nlohmann::json & pt) {

  std::vector<std::string> symbols;
  std::vector<arma::uword> atomic_numbers;
  std::vector<double> coordinates;

  const int charge = pt.at("charge");
  const auto atoms_tree = pt.at("atoms");

  for (const auto & line : atoms_tree) {

    int unit_index = 0;

    if(line.size() != 4) {
      throw Error("incorrect input of geometry detected");
    }

    for(const auto & unit : line) {
      if(unit_index == 0) {

        const auto atom_symbol = unit.get<std::string>();

        symbols.push_back(atom_symbol);

        atomic_numbers.push_back(geometry::periodic_table.at(atom_symbol));

      } else {
        const auto coordinate = unit.get<double>();
        coordinates.push_back(coordinate);
      }

      unit_index++;
    }

  }

  const arma::mat xyz = arma::reshape(arma::vec(coordinates), 3, symbols.size());

  return {symbols, arma::uvec{atomic_numbers}, xyz, charge};

}


}
}