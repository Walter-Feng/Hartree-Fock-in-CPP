#ifndef GEOMETRY_RESOLVE_H
#define GEOMETRY_RESOLVE_H

#include <boost/property_tree/ptree.hpp>
#include "geometry.h"
#include "periodic_table.h"

#include "global/error.h"
#include "util/resolve.h"

namespace hfincpp {
namespace util {

namespace ptree = boost::property_tree;

template<>
geometry::geometry resolve(const ptree::ptree & pt) {

  std::vector<std::string> symbols;
  std::vector<arma::uword> atomic_numbers;
  std::vector<double> coordinates;

  const int charge = pt.get<int>("charge");

  for (const auto & line : pt) {

    int unit_index = 0;

    for(const auto & unit : line.second) {
      if(unit_index == 0) {

        const auto atom_symbol = unit.second.get_value_optional<std::string>();

        if (!atom_symbol) {
          throw Error("Error reading atom symbol");
        }

        symbols.push_back(atom_symbol.value());

        atomic_numbers.push_back(geometry::periodic_table.at(atom_symbol.value()));

      } else {
        const auto coordinate = unit.second.get_value_optional<double>();
        if (!coordinate) {
          throw Error("Error reading atom coordinates");
        }
        coordinates.push_back(coordinate.value());
      }

      unit_index++;
    }

  }

  const arma::mat xyz = arma::reshape(arma::vec(coordinates), 3, symbols.size());

  const geometry::atoms atoms = {symbols, arma::uvec{atomic_numbers}, xyz};

  return {atoms, charge};

}


}
}



#endif //GEOMETRY_RESOLVE_H
