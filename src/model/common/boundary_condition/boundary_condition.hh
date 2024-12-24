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

#ifndef AKANTU_BOUNDARY_CONDITION_HH_
#define AKANTU_BOUNDARY_CONDITION_HH_

#include "aka_common.hh"
#include "boundary_condition_functor.hh"
#include "fe_engine.hh"
#include "mesh.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {

template <class ModelType> class BoundaryCondition {
  /* ------------------------------------------------------------------------ */
  /* Typedefs                                                                 */
  /* ------------------------------------------------------------------------ */
private:
  /* ------------------------------------------------------------------------ */
  /* Constructors / Destructors / Initializers                                */
  /* ------------------------------------------------------------------------ */
public:
  BoundaryCondition() : model(nullptr) {}
  /// Initialize the boundary conditions
  void initBC(ModelType & model, const ID & dof_id, Array<Real> & primal,
              Array<Real> & dual);
  void initBC(ModelType & model, const ID & dof_id, Array<Real> & primal,
              Array<Real> & primal_increment, Array<Real> & dual);
  /* ------------------------------------------------------------------------ */
  /* Methods and accessors                                                    */
  /* ------------------------------------------------------------------------ */
  /// Apply the boundary conditions
  template <typename FunctorType> inline void applyBC(FunctorType && func);

  template <class FunctorType>
  inline void applyBC(FunctorType && func, const std::string & group_name);

  template <class FunctorType>
  inline void applyBC(FunctorType && func, const ElementGroup & element_group);

  AKANTU_GET_MACRO_AUTO_NOT_CONST(Model, *model);
  AKANTU_GET_MACRO_AUTO_NOT_CONST(Primal, *primal);
  AKANTU_GET_MACRO_AUTO_NOT_CONST(Dual, *dual);
  AKANTU_GET_MACRO_AUTO(DOFID, dof_id);

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
  template <class FunctorType,
            BC::Functor::Type type = std::decay_t<FunctorType>::type>
  struct TemplateFunctionWrapper;

private:
  ID dof_id{};
  ModelType * model;
  Array<Real> * primal{nullptr};
  Array<Real> * dual{nullptr};
  Array<Real> * primal_increment{nullptr};
};

} // namespace akantu

#include "boundary_condition_tmpl.hh"

#endif /* AKANTU_BOUNDARY_CONDITION_HH_ */
