#include "catch.h"

#include "integral.h"
//
//TEST_CASE("Check pairwise information") {
//  hfincpp::geometry::Atoms atoms;
//  atoms.atomic_numbers = {8, 1};
//  atoms.xyz = {
//      {0, 0},
//      {0, 1},
//      {0, 0}
//  };
//
//  atoms.symbols = {"O", "H"};
//
//  hfincpp::basis::Basis basis(atoms, "3-21g");
//
//  for(const auto & i : basis.function_labels) {
//    std::cout << i << " ";
//  }
//
//  std::cout << std::endl;
//
//  const auto pair_info =
//      hfincpp::integral::faster_kernel::generate_function_pairs(basis);
//
//  const auto n_pairs = pair_info.exponents.n_elem;
//  CHECK(pair_info.prefactors.n_rows == n_pairs);
//  CHECK(pair_info.centers.n_rows == n_pairs);
//  CHECK(pair_info.screening_conditions.n_rows == n_pairs);
//  CHECK(pair_info.angular_momenta.n_rows == n_pairs);
//  CHECK(pair_info.indexing.n_rows == n_pairs);
//
////  ARMA_DEBUG(pair_info.exponents);
////  ARMA_DEBUG(pair_info.prefactors);
////  ARMA_DEBUG(pair_info.centers);
////  ARMA_DEBUG(pair_info.screening_conditions);
////  ARMA_DEBUG(pair_info.angular_momenta);
////  ARMA_DEBUG(pair_info.indexing);
//}