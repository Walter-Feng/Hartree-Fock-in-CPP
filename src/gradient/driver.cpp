#include "driver.h"

#include "hf/rhf.h"

namespace hfincpp::gradient {

GradientDriver driver(const nlohmann::json & input,
                      const geometry::Atoms & atoms,
                      const basis::Basis & basis,
                      const std::string method) {
  if(method == "rhf") {
    return hf::gradient::driver(input, atoms, basis);
  } else {
    throw Error("method not recognized");
  }
}

}