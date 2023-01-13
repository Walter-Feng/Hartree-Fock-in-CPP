#include "run.h"
#include "util/json.h"
#include "basis/basis.h"
#include "geometry/resolve.h"
#include "gradient/gradient.h"
#include "hf/rhf.h"

namespace hfincpp {

nlohmann::json run(const nlohmann::json & input) {

  const geometry::Atoms atoms =
      util::resolve<geometry::Atoms>(input.at("geometry"));
  const std::string basis_string = input.at("basis");
  const basis::Basis basis(atoms, basis_string);
  const std::string method = input.at("method");


  nlohmann::json output;
  output["input"] = input;

  if(method == "rhf") {
    output["output"] = hf::rhf(input, atoms, basis);
  }

  if(input.contains("extra_commands")) {
    const nlohmann::json extra_commands = input.at("extra_commands");

    for(const auto & command : extra_commands) {
      if(command == "gradient") {
        output["gradient"] =
            gradient::gradient(input, atoms, basis, method)["gradient"];
      }
    }
  }

  return output;
}

}
