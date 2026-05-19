#===============================================================================
# Copyright (©) 2011-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
# Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
#
# This file is part of Akantu
# 
# Akantu is free software: you can redistribute it and/or modify it under the
# terms of the GNU Lesser General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option) any
# later version.
# 
# Akantu is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
# details.
# 
# You should have received a copy of the GNU Lesser General Public License along
# with Akantu. If not, see <http://www.gnu.org/licenses/>.
#
#===============================================================================


package_declare(structural_mechanics DEFAULT OFF
  DESCRIPTION "Use Structural mechanics model package of Akantu"
  DEPENDS implicit)

package_declare_sources(structural_mechanics
  fe_engine/element_class_structural.hh
  fe_engine/element_classes/element_class_bernoulli_beam_inline_impl.hh
  fe_engine/element_classes/element_class_kirchhoff_shell_inline_impl.hh
  fe_engine/element_classes/element_class_hermite_inline_impl.hh
  fe_engine/fe_engine_template_tmpl_struct.hh
  fe_engine/shape_structural.cc
  fe_engine/shape_structural.hh
  fe_engine/shape_structural_inline_impl.hh
  io/mesh_io/mesh_io_msh_struct.cc
  io/mesh_io/mesh_io_msh_struct.hh
  model/structural_mechanics/structural_elements/structural_element_bernoulli_beam_2.hh
  model/structural_mechanics/structural_elements/structural_element_bernoulli_beam_3.hh
  model/structural_mechanics/structural_elements/structural_element_kirchhoff_shell.hh
  model/structural_mechanics/structural_mechanics_model.cc
  model/structural_mechanics/structural_mechanics_model.hh
  model/structural_mechanics/structural_mechanics_model_boundary.cc
  model/structural_mechanics/structural_mechanics_model_inline_impl.hh
  model/structural_mechanics/structural_mechanics_model_mass.cc
  )

