#ifndef EXE_UTIL_PTREE_H
#define EXE_UTIL_PTREE_H

#include <armadillo>
#include <boost/property_tree/ptree.hpp>
#include "global/error.h"

namespace hfincpp {
namespace util {

namespace ptree = boost::property_tree;

template<typename T>
struct std_string_id_translater {
  typedef T internal_type;
  typedef T external_type;

  boost::optional<T> get_value(const T & v) {
    return v.substr(1, v.size() - 2);
  }

  boost::optional<T> put_value(const T & v) { return '"' + v + '"'; }
};

template<typename T>
std::vector<T> get_list(const ptree::ptree & pt) {

  std::vector<T> result;

  for (const auto & unit : pt) {

    const auto value = unit.second.get_value_optional<T>();
    if (!value) {
      throw Error("Error reading value");
    }
    result.push_back(unit.second.get_value<T>());
  }

  return result;
}

template<typename Output_type>
std::vector<Output_type> get_list(const ptree::ptree & pt,
                                  const std::function<Output_type(const ptree::ptree &)> & parser) {

  std::vector<Output_type> result;

  for (const auto & unit : pt) {

    const std::optional<Output_type> value = parser(unit.second);
    if (!value) {
      throw Error("Error reading value");
    }
    result.push_back(value.value());
  }

  return result;
}

template<typename T>
arma::Mat<T> get_mat(const ptree::ptree & list_tree) {

  arma::Mat<T> result;

  for (const auto & line : list_tree) {
    const auto list_elements = arma::Col<T>(get_list<T>(line.second));
    result = arma::join_rows(result, list_elements);
  }

  return result;

}

template<typename Output_type>
arma::field<Output_type> get_mat_object(const ptree::ptree & list_tree,
                                        const std::function<Output_type(const ptree::ptree &)> & parser) {

  std::vector<std::vector<Output_type>> result_in_std_vector;

  for (const auto line : list_tree) {
    const auto list_elements = get_list(line.second, parser);
    result_in_std_vector.push_back(list_elements);
  }

  arma::field<Output_type> result(result_in_std_vector.size(), result_in_std_vector[0].size());

#pragma omp parallel for
  for(arma::uword i=0; i<result.n_rows; i++) {
    for(arma::uword j=0; j<result.n_cols; j++) {
      result(i,j) = result_in_std_vector[i][j];
    }
  }

  return result;

}


template<typename T>
arma::Cube<T> get_cube(const ptree::ptree & list_tree) {

  arma::Cube<T> result;

  for (const auto & line1 : list_tree) {

    const auto list_of_vecs = get_mat<T>(line1.second);
    arma::join_slices(result, list_of_vecs);
  }

  return result;
}

template<typename T>
auto put(ptree::ptree & result, const std::string & path, const T & value) {
  return result.put<T>(path, value);
}

template<>
inline
auto
put(ptree::ptree & result, const std::string & path, const arma::vec & value) {

  ptree::ptree head;

  const auto write = [&head](const double & element) {

    ptree::ptree child;

    child.put<double>("", element);
    head.push_back(std::make_pair("", child));
  };

  value.for_each(write);

  return result.put_child(path, head);
}

template<>
inline
auto
put(ptree::ptree & result, const std::string & path, const arma::mat & value) {

  ptree::ptree head;

  const auto write = [&head](const arma::vec & col) {

    ptree::ptree child_for_col;

    const auto write_in_col = [&child_for_col](const double & element) {
      ptree::ptree child_in_col;
      child_in_col.put<double>("", element);
      child_for_col.push_back(std::make_pair("", child_in_col));
    };

    col.for_each(write_in_col);

    head.push_back(std::make_pair("", child_for_col));
  };

  value.each_col(write);

  return result.put_child(path, head);
}

template<>
inline
auto
put(ptree::ptree & result, const std::string & path, const arma::cube & value) {

  ptree::ptree head;

  const auto write = [&head](const arma::mat & element) {

    ptree::ptree child;

    put(child, "", element);
    head.push_back(std::make_pair("", child));
  };

  value.each_slice(write);

  return result.put_child(path, head);
}

template<>
inline
auto put(ptree::ptree & result, const std::string & path,
         const arma::cx_vec & value) {

  ptree::ptree head;

  const arma::vec real = arma::real(value);
  const arma::vec imag = arma::imag(value);
  put(head, "real", real);
  put(head, "imag", imag);

  return result.put_child(path, head);

}

template<>
inline
auto put(ptree::ptree & result, const std::string & path,
         const arma::cx_mat & value) {

  ptree::ptree head;

  const arma::vec real = arma::real(value);
  const arma::vec imag = arma::imag(value);

  put(head, "real", real);
  put(head, "imag", imag);

  return result.put_child(path, head);
}

template<>
inline
auto put(ptree::ptree & result, const std::string & path,
         const arma::cx_cube & value) {

  ptree::ptree head;

  put(head, "real", arma::real(value));
  put(head, "imag", arma::imag(value));

  return result.put_child(path, head);
}

template<>
inline
auto put(ptree::ptree & result, const std::string & path,
         const std::string & value) {

  return
      result.put<std::string>(path, value,
                              std_string_id_translater<std::string>());
}

template<typename T>
inline
auto put(ptree::ptree & result, const std::string & path,
         const std::vector<T> & value) {
  ptree::ptree head;

  const auto write = [&head](const T element) {

    ptree::ptree child;

    child.put<T>("", element);
    head.push_back(std::make_pair("", child));
  };

  for (const auto element : value) {
    write(element);
  }

  return result.put_child(path, head);
}

template<>
inline
auto put(ptree::ptree & result, const std::string & path,
         const std::vector<std::string> & value) {
  ptree::ptree head;

  const auto write = [&head](const std::string & element) {

    ptree::ptree child;

    put(child, "", element);
    head.push_back(std::make_pair("", child));
  };

  for (const auto & element : value) {
    write(element);
  }

  return result.put_child(path, head);
}


}
}

#endif //EXE_UTIL_PTREE_H
