#ifndef UTIL_PRINTER_H
#define UTIL_PRINTER_H

#include <fmt/format.h>

#include "cx_double.h"

namespace hfincpp {
template<typename T>
std::string format(const arma::Mat<T> & arma,
                   const int precision = 6,
                   const int width = 10,
                   const std::string aligned = ">") {

  std::string result = "";

  for (arma::uword i = 0; i < arma.n_rows; i++) {
    for (arma::uword j = 0; j < arma.n_cols; j++) {
      const auto formatted =
          fmt::format("{:.{}g}", arma(i, j), precision);
      result += fmt::format("{:" + aligned + "{}}", formatted, width);
    }
    result += fmt::format("\n");
  }

  return result;
}

template<typename T>
std::string format(const T number,
                   const int precision = 6,
                   const int width = -1,
                   const std::string aligned = ">") {

  if(width < 0) {
    return fmt::format("{}",number);
  } else {
    const auto formatted =
        fmt::format("{:.{}}", number, precision);

    return fmt::format("{:" + aligned + "{}}", formatted, width);
  }
}

template<>
inline
std::string format(const double number,
                   const int precision,
                   const int width,
                   const std::string aligned) {
  const auto formatted =
      fmt::format("{:.{}g}", number, precision);

  if(width < 0) {
    return formatted;
  } else {
    return fmt::format("{:" + aligned + "{}}", formatted, width);
  }
}

template<>
inline
std::string format(const cx_double number,
                   const int precision,
                   const int width,
                   const std::string aligned) {
  const auto formatted =
      "(" + fmt::format("{:.{}g}", std::real(number), precision) + "," +
      fmt::format("{:.{}g}", std::imag(number), precision) + ")";

  if(width < 0) {
    return formatted;
  } else {
    return fmt::format("{:" + aligned + "{}}", formatted, width);
  }
}


template<typename T>
void print(const arma::Row<T> & arma,
           const int width = 10,
           const int precision = 6,
           const std::string aligned = ">") {
  for (arma::uword j = 0; j < arma.n_cols; j++) {
    const std::string formatted = format<T>(arma(j), precision);
    fmt::print("{:" + aligned + "{}}", formatted, width);
  }
}

template<typename T>
void print(const arma::Col<T> & arma,
           const int width = 10,
           const int precision = 6,
           const std::string aligned = ">") {
  for (arma::uword i = 0; i < arma.n_rows; i++) {
    const std::string formatted = format<T>(arma(i), precision);
    fmt::print("{:" + aligned + "{}}", formatted, width);
    fmt::print("\n");
  }
}

template<typename T>
void print(const arma::Mat<T> & arma,
           const int width = 10,
           const int precision = 6,
           const std::string aligned = ">") {
  for (arma::uword i = 0; i < arma.n_rows; i++) {
    for (arma::uword j = 0; j < arma.n_cols; j++) {
      const std::string formatted =
          fmt::format("{:.{}}", arma(i, j), precision);
      fmt::print("{:" + aligned + "{}}", formatted, width);
    }
    fmt::print("\n");
  }
}

template<typename T>
void print(const T number,
           const int width = 10,
           const int precision = 6,
           const std::string aligned = ">") {
  const auto formatted = format<T>(number, precision);
  fmt::print("{:" + aligned + "{}}", formatted, width);
}


template<typename State>
using Printer = std::function<int(const State & state,
                                  const double computation_time,
                                  const int iter,
                                  const int print_level,
                                  const bool print_header)>;

template<typename State>
Printer<State> mute = [](const State &,
                         const arma::uword,
                         const double,
                         const int,
                         const bool) -> int {

  return 0;
};

template<typename State>
Printer<State> operator<<(const Printer<State> a, const Printer<State> b) {
  return [a, b](const State & state,
                const arma::uword index,
                const double time,
                const int print_level = 1,
                const bool print_header = false) -> int {

    const int a_length = a(state, index, time, print_level, print_header);
    const int b_length = b(state, index, time, print_level, print_header);

    return std::max(a_length, b_length);

  };
}

}

#endif //UTIL_PRINTER_H
