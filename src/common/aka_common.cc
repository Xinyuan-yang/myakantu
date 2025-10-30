/**
 * Copyright (©) 2010-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * This file is part of Akantu
 *
 * Akantu is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * Akantu is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
 * details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with Akantu. If not, see <http://www.gnu.org/licenses/>.
 */

/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "aka_fwd.hh"
#include "aka_random_generator.hh"
#include "communicator.hh"

#include "cppargparse.hh"
#include "parser.hh"

#include "communication_tag.hh"
/* -------------------------------------------------------------------------- */
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <iomanip>
#include <ios>
#include <limits>
#include <memory>
#include <ostream>
#include <string>
#include <type_traits>
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
void initialize(int & argc, char **& argv) {
  AKANTU_DEBUG_IN();

  initialize("", argc, argv);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void initialize(const std::string & input_file) {
  int argc = 0;
  char ** argv = nullptr;

  initialize(input_file, argc, argv);
}

/* -------------------------------------------------------------------------- */
void initialize(const std::string & input_file, int & argc, char **& argv) {
  AKANTU_DEBUG_IN();
  Communicator & comm = Communicator::getStaticCommunicator(argc, argv);

  Tag::setMaxTag(comm.getMaxTag());

  debug::debugger.setParallelContext(comm.whoAmI(), comm.getNbProc());
  debug::setDebugLevel(dblError);

  static_argparser.setParallelContext(comm.whoAmI(), comm.getNbProc());
  static_argparser.setExternalExitFunction(debug::exit);
  static_argparser.addArgument("--aka_input_file", "Akantu's input file", 1,
                               cppargparse::_string, std::string());
  static_argparser.addArgument(
      "--aka_debug_level",
      std::string("Akantu's overall debug level") +
          std::string(" (0: error, 1: exceptions, 4: warnings, 5: info, ..., "
                      "100: dump") +
          std::string(" more info on levels can be found in aka_error.hh)"),
      1, cppargparse::_integer, (long int)(dblWarning));

  static_argparser.addArgument(
      "--aka_print_backtrace",
      "Should Akantu print a backtrace in case of error", 0,
      cppargparse::_boolean, false, true);

  static_argparser.addArgument("--aka_seed", "The seed to use on prank 0", 1,
                               cppargparse::_integer);

  static_argparser.parse(argc, argv, cppargparse::_remove_parsed);

  std::string infile = static_argparser["aka_input_file"];
  if (infile.empty()) {
    infile = input_file;
  }
  debug::debugger.printBacktrace(static_argparser["aka_print_backtrace"]);

  if (not infile.empty()) {
    readInputFile(infile);
  }

  long int seed;
  char * env_seed = std::getenv("AKA_SEED");
  if (env_seed != nullptr) {
    seed = std::atol(env_seed);
  } else if (static_argparser.has("aka_seed")) {
    seed = static_argparser["aka_seed"];
  } else {
    seed =
        static_parser.getParameter("seed", time(nullptr), _ppsc_current_scope);
  }

  seed *= (comm.whoAmI() + 1);
  RandomGenerator<Idx>::seed(seed);

  long int dbl_level = static_argparser["aka_debug_level"];
  debug::setDebugLevel(DebugLevel(dbl_level));

  char * env_debug_level = std::getenv("AKA_DEBUG_LEVEL");
  if (env_debug_level != nullptr) {
    debug::setDebugLevel(DebugLevel(std::atoi(env_debug_level)));
  }

  AKANTU_DEBUG_INFO("Random seed set to " << seed);

  std::atexit(finalize);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void finalize() {}

/* -------------------------------------------------------------------------- */
void readInputFile(const std::string & input_file) {
  static_parser.parse(input_file);
}

/* -------------------------------------------------------------------------- */
cppargparse::ArgumentParser & getStaticArgumentParser() {
  return static_argparser;
}

/* -------------------------------------------------------------------------- */
Parser & getStaticParser() { return static_parser; }

/* -------------------------------------------------------------------------- */
const ParserSection & getUserParser() {
  return *(static_parser.getSubSections(ParserType::_user).first);
}

std::unique_ptr<Communicator> Communicator::static_communicator;

std::ostream & operator<<(std::ostream & stream, NodeFlag flag) {
  using under = std::underlying_type_t<NodeFlag>;
  auto digits = static_cast<int>(
      std::log(std::numeric_limits<under>::max() + 1) / std::log(16));
  std::ios_base::fmtflags ff;
  ff = stream.flags();
  auto value = static_cast<std::common_type_t<under, unsigned int>>(flag);
  stream << "0x" << std::hex << std::setw(digits) << std::setfill('0') << value;
  stream.flags(ff);
  return stream;
}

} // namespace akantu
