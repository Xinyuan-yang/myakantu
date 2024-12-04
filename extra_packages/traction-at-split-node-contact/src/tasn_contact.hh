/**
 * Copyright (©) 2016-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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

// ast common
#include "manual_restart.hh"
#include "parameter_reader.hh"
#include "synchronized_array.hh"

// functions
#include "boundary_functions.hh"
#include "node_filter.hh"

// boundary conditions
#include "force_based_dirichlet.hh"
#include "inclined_flat_dirichlet.hh"
#include "spring_bc.hh"

// ntn/ntrf contact
#include "mIIasym_contact.hh"
#include "ntn_base_contact.hh"
#include "ntn_contact.hh"
#include "ntrf_contact.hh"

// ntn/ntrf friction
#include "ntn_base_friction.hh"
#include "ntn_friction.hh"
#include "ntrf_friction.hh"

// friction regularisations
#include "ntn_fricreg_no_regularisation.hh"
#include "ntn_fricreg_rubin_ampuero.hh"
#include "ntn_fricreg_simplified_prakash_clifton.hh"

// friction laws
#include "ntn_friclaw_coulomb.hh"
#include "ntn_friclaw_linear_cohesive.hh"
#include "ntn_friclaw_linear_slip_weakening.hh"
#include "ntn_friclaw_linear_slip_weakening_no_healing.hh"
#include "ntn_friclaw_rate_and_state.hh"
#include "ntn_friclaw_rate_and_state_dieterich.hh"
#include "ntn_friclaw_rate_and_state_flashheating.hh"
#include "ntn_friclaw_rate_and_state_vw.hh"

// initiation of friction
#include "ntn_initiation_function.hh"
