#include "global/global.h"

#include <args.hxx>
#include <armadillo>
#include <iostream>
#include <fmt/core.h>
#include <json.hpp>

#include "run.h"
#include "global/error.h"
#include "util/time.h"


int main(const int argc, const char * argv[]) {

  using namespace hfincpp;

  Timer global_time;

  args::ArgumentParser parser("This is the executable of a Hartree-Fock Program"
                              "using full C++ standard."
                              "Written by Walter Feng (github.com/Walter-Feng)");

  args::HelpFlag help(parser, "help",
                      "Display this help menu", {'h', "help"});

  args::Positional<std::string> input_flag(parser, "input",
                                           "The input file (in json format)");

  args::ValueFlag<std::string> str_input(parser, "string",
                                         "String form of json input", {'s'});

  args::Flag print_json(parser, "json",
                        "Print as json format."
                        " Remember to set print_level = -1 if you would like to "
                        "do hfincpp input.json > output.json", {'j'});


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



  nlohmann::json input;
  if(str_input) {
    input = nlohmann::json::parse(args::get(str_input));
  } else {
    const std::string input_filename = args::get(input_flag);
    if(input_filename.empty()) {
      throw Error("No string input or file input is given");
    }
    std::ifstream input_file_stream(args::get(input_flag));

    input_file_stream >> input;
  }

  nlohmann::json result = run(input);

  if(input["print_level"] >= 1) {
    fmt::print("Total time elapsed: {} s\n", global_time.elapsed());
  }

  if(print_json) {
    std::cout << std::setw(4) << result << std::endl;
  }


  return 0;

}