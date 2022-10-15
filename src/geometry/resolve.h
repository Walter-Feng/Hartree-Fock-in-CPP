#ifndef GEOMETRY_RESOLVE_H
#define GEOMETRY_RESOLVE_H

#include <json.hpp>
#include "geometry.h"
#include "periodic_table.h"

#include "global/error.h"
#include "util/resolve.h"

namespace hfincpp {
namespace util {

template<>
geometry::Geometry resolve(const nlohmann::json & pt) {

  std::vector<std::string> symbols;
  std::vector<arma::uword> atomic_numbers;
  std::vector<double> coordinates;

  const int charge = pt.at("charge");

  for (const auto & line : pt) {

    int unit_index = 0;

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

  const geometry::Atoms atoms = {symbols, arma::uvec{atomic_numbers}, xyz};

  return {atoms, charge};

}


}
}



#endif //GEOMETRY_RESOLVE_H
