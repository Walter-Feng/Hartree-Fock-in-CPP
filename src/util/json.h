#ifndef UTIL_JSON_H
#define UTIL_JSON_H

#include <armadillo>
#include <nlohmann/json.hpp>
#include "global/error.h"

namespace hfincpp {
namespace util {

template<typename T>
std::vector<T> get_list(const nlohmann::json & pt) {

  std::vector<T> result;

  for (const auto & unit : pt) {
    result.push_back(unit.value<T>());
  }

  return result;
}

template<typename T>
arma::Mat<T> get_mat(const nlohmann::json & list_tree) {

  arma::Mat<T> result;

  for (const auto & line : list_tree) {
    const auto list_elements = arma::Col<T>(get_list<T>(line));
    result = arma::join_rows(result, list_elements);
  }

  return result;

}

template<typename T>
auto put(nlohmann::json & result, const std::string & path, const T & value) {
  (result[path] = value);
}

template<>
inline
auto
put(nlohmann::json & result, const std::string & path, const arma::vec & value) {

  std::vector<double> converted_array;

  for(arma::uword i=0; i<value.n_elem; i++) {
    converted_array.push_back(value(i));
  }

  result[path] = converted_array;
}

template<>
inline
auto
put(nlohmann::json & result, const std::string & path, const arma::mat & value) {

  std::vector<std::vector<double>> converted_matrix;

  for(arma::uword i=0; i<value.n_rows; i++) {
    std::vector<double> row_vector;
    for(arma::uword j=0; j<value.n_cols; j++) {
      row_vector.push_back(value(i, j));
    }
    converted_matrix.push_back(row_vector);
  }

  (result[path] = converted_matrix);
}

template<>
inline
auto put(nlohmann::json & result, const std::string & path,
         const arma::cx_vec & value) {

  nlohmann::json head;
  const arma::vec real = arma::real(value);
  const arma::vec imag = arma::imag(value);
  put(head, "real", real);
  put(head, "imag", imag);

  result[path] = head;
}

template<>
inline
auto put(nlohmann::json & result, const std::string & path,
         const arma::cx_mat & value) {

  nlohmann::json head;

  const arma::vec real = arma::real(value);
  const arma::vec imag = arma::imag(value);

  put(head, "real", real);
  put(head, "imag", imag);

  result[path] = head;
}


template<typename T>
inline
auto put(nlohmann::json & result, const std::string & path,
         const std::vector<T> & value) {
  result[path] = value;
}

}
}

#endif //UTIL_JSON_H
