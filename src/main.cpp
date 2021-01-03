#include "global/global.h"

#include <args.hxx>
#include <armadillo>
#include <iostream>
#include <fmt/core.h>
#include <nlohmann/json.hpp>

#include "run.h"
#include "util/time.h"


int main(const int argc, const char * argv[]) {

  using namespace hfincpp;

  Timer global_time;

  args::ArgumentParser parser("This is the executable of a Hartree-Fock Program"
                              "using full C++ standard."
                              "Written by Walter Feng (github.com/Walter-Feng)");

  args::HelpFlag help(parser, "help", "Display this help menu", {'h', "help"});

  args::Positional<std::string> input_flag(parser, "input",
                                           "The input file (in json format)");

  try {
    parser.ParseCLI(argc, argv);
  }
  catch (const args::Help &) {
    std::cout << parser << std::endl;
    return 0;
  }
  catch (const args::ParseError & e) {
    std::cout << e.what() << std::endl;
    std::cout << parser << std::endl;
    return 1;
  }
  catch (const args::ValidationError & e) {
    std::cout << e.what() << std::endl;
    std::cout << parser << std::endl;
    return 1;
  }

  ///////////////////// Read Input File /////////////////////


  std::ifstream input_file_stream(args::get(input_flag));
  nlohmann::json input;
  input_file_stream >> input;

  nlohmann::json result = run(input);

  fmt::print("Total time elapsed: {} s\n", global_time.elapsed());

  return 0;

}