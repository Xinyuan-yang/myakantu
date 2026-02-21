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
#include "aka_array.hh"
#include "aka_common.hh"
#include "aka_math.hh"
#include "aka_named_argument.hh"
#include "aka_random_generator.hh"
#include "communication_tag.hh"
#include "cppargparse.hh"
#include "parser.hh"
#include "solid_mechanics_model.hh"
#include <string>
#if defined(AKANTU_COHESIVE_ELEMENT)
#include "solid_mechanics_model_cohesive.hh"
#endif
/* -------------------------------------------------------------------------- */
#include <iostream>
#include <limits>

namespace akantu {

/* -------------------------------------------------------------------------- */
/* error.hpp variables                                                        */
/* -------------------------------------------------------------------------- */
namespace debug {
/** \todo write function to get this
 *   values from the environment or a config file
 */
/// standard output for debug messages
std::ostream * _akantu_debug_cout = &std::cerr;

/// standard output for normal messages
std::ostream & _akantu_cout = std::cout;

/// parallel context used in debug messages
std::string _parallel_context;

Debugger debugger;

} // namespace debug

/* -------------------------------------------------------------------------- */
/// Paser for commandline arguments
::cppargparse::ArgumentParser static_argparser;

/// Parser containing the information parsed by the input file given to initFull
Parser static_parser;

bool Parser::permissive_parser = false;

/* -------------------------------------------------------------------------- */
Real Math::tolerance = 1e2 * std::numeric_limits<Real>::epsilon();

/* -------------------------------------------------------------------------- */
const Int _all_dimensions [[gnu::unused]] = Int(-1);

/* -------------------------------------------------------------------------- */
/// Szie of one to be sure to have a non nullptr to compare
const Array<Int> empty_filter(1, 1, "empty_filter");

/* -------------------------------------------------------------------------- */
template <> long int RandomGenerator<Idx>::_seed = 5489;
template <> long int RandomGenerator<UInt>::_seed = 5489U;

template <> std::default_random_engine RandomGenerator<Idx>::generator(5489);
template <> std::default_random_engine RandomGenerator<UInt>::generator(5489U);

/* -------------------------------------------------------------------------- */
int Tag::max_tag = 0;
/* -------------------------------------------------------------------------- */

} // namespace akantu
