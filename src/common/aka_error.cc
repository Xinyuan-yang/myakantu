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
#include "aka_error.hh"
#include "aka_common.hh"
#include "aka_config.hh"
#include "aka_iterators.hh"
#include "aka_random_generator.hh"
/* -------------------------------------------------------------------------- */
#include <boost/preprocessor/stringize.hpp>
#include <cstdio>
#include <cstdlib>
#include <exception>
#include <iostream>
#include <ratio>
#include <sstream>
#include <stdio.h>
#include <string>
#include <vector>

#if (defined(READLINK_COMMAND) || defined(ADDR2LINE_COMMAND)) &&               \
    (!defined(_WIN32))
#include <execinfo.h>
#include <sys/wait.h>
#endif

#include <chrono>
#include <cmath>
#include <cstring>
#include <cxxabi.h>
#include <fstream>
#include <iomanip>
#include <map>
#include <sys/types.h>
#include <unistd.h>

#ifdef AKANTU_USE_MPI
#include <mpi.h>
#endif

/* -------------------------------------------------------------------------- */
namespace akantu {
namespace debug {
/* ------------------------------------------------------------------------ */
std::string demangle(const char * symbol) {
  int status{};
  std::string result;

  if (auto * demangled_name =
          abi::__cxa_demangle(symbol, nullptr, nullptr, &status);
      demangled_name != nullptr) {
    result = demangled_name;
    std::free(demangled_name); // NOLINT
  } else {
    result = symbol;
  }

  return result;
}

/* ------------------------------------------------------------------------ */
// NOLINTBEGIN(cppcoreguidelines-avoid-c-arrays,
// cppcoreguidelines-pro-bounds-pointer-arithmetic,
// performance-inefficient-string-concatenation,
// cppcoreguidelines-owning-memory)
#if (defined(READLINK_COMMAND) || defined(ADDR2LINK_COMMAND)) &&               \
    (!defined(_WIN32))
std::string exec(const std::string & cmd) {
  FILE * pipe = popen(cmd.c_str(), "r");
  if (pipe == nullptr) {
    return "";
  }
  char buffer[1024];
  std::string result;
  while (feof(pipe) == 0) {
    if (fgets(buffer, 128, pipe) != nullptr) {
      result += buffer;
    }
  }

  result = result.substr(0, result.size() - 1);
  pclose(pipe);
  return result;
}
#endif

auto getBacktrace() -> std::vector<std::string> {
  std::vector<std::string> backtrace_lines;
#if not defined(_WIN32)
#if defined(READLINK_COMMAND) && defined(ADDR2LINE_COMMAND)
  std::string me;
  char buf[1024];
  /* The manpage says it won't null terminate.  Let's zero the buffer. */
  memset(buf, 0, sizeof(buf));
  /* Note we use sizeof(buf)-1 since we may need an extra char for NUL. */
  if (readlink("/proc/self/exe", buf, sizeof(buf) - 1) != 0) {
    me = std::string(buf);
  }

  std::ifstream inmaps;
  inmaps.open("/proc/self/maps");
  std::map<std::string, size_t> addr_map;
  std::string line;
  while (inmaps.good()) {
    std::getline(inmaps, line);
    std::stringstream sstr(line);

    size_t first = line.find('-');
    std::stringstream sstra(line.substr(0, first));
    size_t addr{};
    sstra >> std::hex >> addr;

    std::string lib;
    sstr >> lib;
    sstr >> lib;
    sstr >> lib;
    sstr >> lib;
    sstr >> lib;
    sstr >> lib;
    if (not lib.empty() and (addr_map.find(lib) == addr_map.end())) {
      addr_map[lib] = addr;
    }
  }

  if (not me.empty()) {
    addr_map[me] = 0;
  }
#endif

  /// \todo for windows this part could be coded using CaptureStackBackTrace
  /// and SymFromAddr
  const size_t max_depth = 100;
  size_t stack_depth{};
  void * stack_addrs[max_depth];
  char ** stack_strings{nullptr};

  size_t i{};
  stack_depth = backtrace(stack_addrs, max_depth);
  stack_strings = backtrace_symbols(stack_addrs, int(stack_depth));

  /// -1 to remove the call to the printBacktrace function
  for (i = 1; i < stack_depth; i++) {
    std::string bt_line(stack_strings[i]);
    size_t first{};
    size_t second{};

    if ((first = bt_line.find('(')) != std::string::npos &&
        (second = bt_line.find('+')) != std::string::npos) {
      std::string location = bt_line.substr(0, first);
#if defined(READLINK_COMMAND)
      std::string location_cmd =
          std::string(BOOST_PP_STRINGIZE(READLINK_COMMAND)) +
                      std::string(" -f ") + location;
      location = exec(location_cmd);
#endif
      std::string call =
          demangle(bt_line.substr(first + 1, second - first - 1).c_str());
      size_t f = bt_line.find('[');
      size_t s = bt_line.find(']');
      std::string address = bt_line.substr(f + 1, s - f - 1);
      std::stringstream sstra(address);
      size_t addr{};
      sstra >> std::hex >> addr;
      std::string trace = location + " [" + call + "]";
#if defined(READLINK_COMMAND) && defined(ADDR2LINE_COMMAND)
      auto it = addr_map.find(location);
      if (it != addr_map.end()) {
        std::stringstream syscom;
        syscom << BOOST_PP_STRINGIZE(ADDR2LINE_COMMAND)
                                     << " 0x" << std::hex << (addr - it->second)
                                     << " -i -e " << location;
        std::string line = exec(syscom.str());
        trace += " (" + line + ")";
      } else {
#endif
        std::stringstream sstr_addr;
        sstr_addr << std::hex << addr;
        trace += " (0x" + sstr_addr.str() + ")";
#if defined(READLINK_COMMAND) && defined(ADDR2LINE_COMMAND)
      }
#endif
      backtrace_lines.push_back(trace);
    } else {
      backtrace_lines.push_back(bt_line);
    }
  }

  free(stack_strings);

#endif
  return backtrace_lines;
}
// NOLINTEND(cppcoreguidelines-avoid-c-arrays,
// cppcoreguidelines-pro-bounds-pointer-arithmetic,
// performance-inefficient-string-concatenation,
// cppcoreguidelines-owning-memory)
/* ------------------------------------------------------------------------ */
void printBacktrace(const std::vector<std::string> & backtrace) {
  auto w = size_t(std::floor(std::log10(double(backtrace.size()))) + 1);
  std::cerr << "BACKTRACE :  " << backtrace.size() << " stack frames.\n";
  for (auto && data : enumerate(backtrace)) {
    std::cerr << "  [" << std::setw(int(w)) << (std::get<0>(data) + 1) << "] "
              << std::get<1>(data) << "\n";
  }
  std::cerr << "END BACKTRACE"
            << "\n";
}

/* ------------------------------------------------------------------------ */
namespace {
void terminate_handler() {
  auto eptr = std::current_exception();
  auto * t = abi::__cxa_current_exception_type();
  auto name = (t != nullptr) ? demangle(t->name()) : std::string("unknown");
  try {
    if (eptr) {
      std::rethrow_exception(eptr);
    } else {
      printBacktrace();
      std::cerr << AKANTU_LOCATION
                << "!! Execution terminated for unknown reasons !!"
                << "\n";
    }
  } catch (Exception & e) {
    printBacktrace(e.backtrace());
    std::cerr << "!! Uncaught akantu::Exception of type " << name
              << " !!\nwhat(): \"" << e.what() << "\""
              << "\n";
  } catch (std::exception & e) {
    std::cerr << "!! Uncaught exception of type " << name << " !!\nwhat(): \""
              << e.what() << "\""
              << "\n";
  } catch (...) {
    std::cerr << "!! Something strange of type \"" << name
              << "\" was thrown.... !!"
              << "\n";
  }

  if (debugger.printBacktrace()) {
    std::cerr << "Random generator seed: " << RandomGenerator<Int>::seed()
              << "\n";
    printBacktrace();
  }
}
} // namespace

/* ------------------------------------------------------------------------ */
/* ------------------------------------------------------------------------ */
Debugger::Debugger() noexcept : cout(&std::cerr) {
  // initSignalHandler();
  std::set_terminate(terminate_handler);
}

/* ------------------------------------------------------------------------ */
Debugger::~Debugger() {
  if (file_open) {
    dynamic_cast<std::ofstream *>(cout)->close();
    delete cout;
  }
}

/* ------------------------------------------------------------------------ */
void Debugger::exit(int status) {
  if (status != 0) {
    std::terminate();
  }

  std::exit(0);
}

/*------------------------------------------------------------------------- */
void Debugger::throwException(const std::string & info,
                              const std::string & file, unsigned int line,
                              __attribute__((unused)) bool silent,
                              __attribute__((unused))
                              const std::string & location,
                              const std::string & module) const
    noexcept(false) {

#if !defined(AKANTU_NDEBUG)
  if (not silent) {
    printMessage("###", dblWarning, info + " " + location, module);
  }
#endif

  debug::Exception ex(info, file, line);
  ex.setModule(module);
  throw ex;
}

/* ------------------------------------------------------------------------ */
void Debugger::printMessage(const std::string & prefix,
                            const DebugLevel & level, const std::string & info,
                            const std::string & module) const {
  if AKANTU_UNLIKELY (testLevel(level, module)) {
    double timestamp =
        std::chrono::duration_cast<std::chrono::duration<double, std::micro>>(
            std::chrono::system_clock::now().time_since_epoch())
            .count();

    *(cout) << parallel_context << "{" << (size_t)timestamp << "} " << prefix
            << " " << info << "\n";
  }
}

/* ------------------------------------------------------------------------ */
void Debugger::setDebugLevel(const DebugLevel & level) { this->level = level; }

/* ------------------------------------------------------------------------ */
auto Debugger::getDebugLevel() const -> const DebugLevel & {
  return this->level;
}

/* ------------------------------------------------------------------------ */
void Debugger::setLogFile(const std::string & filename) {
  if (file_open) {
    dynamic_cast<std::ofstream *>(cout)->close();
    delete cout;
  }

  // NOLINTNEXTLINE(cppcoreguidelines-owning-memory)
  auto * fileout = new std::ofstream(filename.c_str());
  file_open = true;
  cout = fileout;
}

auto Debugger::getOutputStream() -> std::ostream & { return *cout; }

/* ------------------------------------------------------------------------ */
void Debugger::setParallelContext(int rank, int size) {
  std::stringstream sstr;
  Int pad = std::ceil(std::log10(size));
  sstr << "<" << getpid() << ">[R" << std::setfill(' ') << std::right
       << std::setw(pad) << rank << "|S" << size << "] ";
  parallel_context = sstr.str();
}

void setDebugLevel(const DebugLevel & level) { debugger.setDebugLevel(level); }

auto getDebugLevel() -> const DebugLevel & { return debugger.getDebugLevel(); }

/* ------------------------------------------------------------------------ */
void exit(int status) { Debugger::exit(status); }

} // namespace debug
} // namespace akantu
