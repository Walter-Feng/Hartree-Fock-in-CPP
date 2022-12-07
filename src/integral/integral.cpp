#include "integral.h"

namespace hfincpp::integral {

std::vector<GaussianFunction> GaussianFunction::gradient(const int xyz_index) const {

  arma::Col<int>::fixed<3> angular_change =
      arma::Col<int>::fixed<3>(arma::fill::zeros);
  angular_change(xyz_index) = 1;

  std::vector<GaussianFunction> result = {
      {center, angular + angular_change, exponent, -coef * 2.0 * exponent}
  };

  if(angular(xyz_index) > 0) {
    result.push_back(
        {center, angular - angular_change, exponent, coef * angular(xyz_index)}
        );
  }

  return result;
}

std::vector<GaussianFunction> GaussianFunction::laplace() const {
  std::vector<GaussianFunction> result;
  for(int xyz_index=0; xyz_index<3; xyz_index++) {
    const auto first_order = this->gradient(xyz_index);
    for(const auto & i_gto : first_order) {
      const auto second_order = i_gto.gradient(xyz_index);
      result.insert(result.end(), second_order.begin(), second_order.end());
    }
  }

  return result;
}
}
