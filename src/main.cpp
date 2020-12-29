#include "global/global.h"

#include <args.hxx>
#include <armadillo>
#include <iostream>
#include <fmt/core.h>
#include <boost/property_tree/ptree.hpp>
#include "util/json_parser.h"

#include "run.h"
#include "util/time.h"


int main(const int argc, const char * argv[]) {

  using namespace hfincpp;

  Timer global_time;

  namespace ptree = boost::property_tree;

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

  ptree::ptree input;

  ptree::read_json(args::get(input_flag), input);

  ptree::ptree result = run(input);

  if (input.get_optional<std::string>("json")) {
    result.put<double>("time_elapsed", global_time.elapsed());
    ptree::write_json(input.get<std::string>("json"), result);
  }

  fmt::print("Total time elapsed: {} s\n", global_time.elapsed());

  return 0;

}