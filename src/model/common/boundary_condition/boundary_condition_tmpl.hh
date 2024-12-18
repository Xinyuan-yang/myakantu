/**
 * Copyright (©) 2013-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "boundary_condition.hh"
#include "element_group.hh"
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_BOUNDARY_CONDITION_TMPL_HH_
#define AKANTU_BOUNDARY_CONDITION_TMPL_HH_

namespace akantu {

/* -------------------------------------------------------------------------- */
template <typename ModelType>
void BoundaryCondition<ModelType>::initBC(ModelType & model,
                                          Array<Real> & primal,
                                          Array<Real> & dual) {
  this->model = &model;
  this->primal = &primal;
  this->dual = &dual;
}

/* -------------------------------------------------------------------------- */
template <typename ModelType>
void BoundaryCondition<ModelType>::initBC(ModelType & model,
                                          Array<Real> & primal,
                                          Array<Real> & primal_increment,
                                          Array<Real> & dual) {
  this->initBC(model, primal, dual);
  this->primal_increment = &primal_increment;
}

/* -------------------------------------------------------------------------- */
/* Partial specialization for DIRICHLET functors */
template <typename ModelType>
template <typename FunctorType>
struct BoundaryCondition<ModelType>::TemplateFunctionWrapper<
    FunctorType, BC::Functor::_dirichlet> {
  static inline void applyBC(FunctorType && func, const ElementGroup & group,
                             BoundaryCondition<ModelType> & bc_instance) {
    auto & model = bc_instance.getModel();
    auto & primal = bc_instance.getPrimal();
    auto & boundary_flags = model.getBlockedDOFs();

    FunctorType func_(func);

    const auto & coords = model.getMesh().getNodes();
    Int dim = model.getMesh().getSpatialDimension();

    auto it = zip(make_view(primal, primal.getNbComponent()),
                  make_view(boundary_flags, boundary_flags.getNbComponent()),
                  make_const_view(coords, dim))
                  .begin();
    for (auto && n : group.getNodeGroup()) {
      auto && [primal_, flags_, coords_] = it[n];
      // The copy it to avoid the user to template is functor
      auto && primal = primal_;
      auto && flags = flags_;
      auto && coords = coords_;
      func_(n, flags, primal, coords);
    }
  }
};

/* -------------------------------------------------------------------------- */
/* Partial specialization for NEUMANN functors */
template <typename ModelType>
template <typename FunctorType>
struct BoundaryCondition<ModelType>::TemplateFunctionWrapper<
    FunctorType, BC::Functor::_neumann> {
  static inline void applyBC(FunctorType && func, const ElementGroup & group,
                             BoundaryCondition<ModelType> & bc_instance) {
    auto dim = bc_instance.getModel().getSpatialDimension();

    if (dim == 1) {
      AKANTU_TO_IMPLEMENT();
    }
    applyBC(std::forward<decltype(func)>(func), group, bc_instance, _not_ghost);
    applyBC(std::forward<decltype(func)>(func), group, bc_instance, _ghost);
  }

  static inline void applyBC(FunctorType && func, const ElementGroup & group,
                             BoundaryCondition<ModelType> & bc_instance,
                             GhostType ghost_type) {
    auto & model = bc_instance.getModel();
    auto & dual = bc_instance.getDual();
    const auto & mesh = model.getMesh();
    const auto & nodes_coords = mesh.getNodes();
    const auto & fem_boundary = model.getFEEngineBoundary();

    FunctorType func_(func);

    Int dim = model.getSpatialDimension();
    Int nb_degree_of_freedom = dual.getNbComponent();

    IntegrationPoint quad_point;
    quad_point.ghost_type = ghost_type;

    // Loop over the boundary element types
    for (auto && type : group.elementTypes(dim - 1, ghost_type)) {
      const auto & element_ids = group.getElements(type, ghost_type);

      auto nb_quad_points =
          fem_boundary.getNbIntegrationPoints(type, ghost_type);
      auto nb_elements = element_ids.size();
      auto nb_nodes_per_element = mesh.getNbNodesPerElement(type);

      Array<Real> dual_before_integ(nb_elements * nb_quad_points,
                                    nb_degree_of_freedom, 0.);
      Array<Real> quad_coords(nb_elements * nb_quad_points, dim);

      const auto & normals_on_quad =
          fem_boundary.getNormalsOnIntegrationPoints(type, ghost_type);

      fem_boundary.interpolateOnIntegrationPoints(
          nodes_coords, quad_coords, dim, type, ghost_type, element_ids);
      auto normals_begin = normals_on_quad.cbegin(dim);
      decltype(normals_begin) normals_iter;
      auto quad_coords_iter = quad_coords.cbegin(dim);
      auto dual_iter = dual_before_integ.begin(nb_degree_of_freedom);

      quad_point.type = type;
      for (auto el : element_ids) {
        quad_point.element = el;
        normals_iter = normals_begin + el * nb_quad_points;
        for (auto q : arange(nb_quad_points)) {
          quad_point.num_point = q;

          func_(quad_point, *dual_iter, *quad_coords_iter, *normals_iter);

          ++dual_iter;
          ++quad_coords_iter;
          ++normals_iter;
        }
      }

      Array<Real> dual_by_shapes(nb_elements * nb_quad_points,
                                 nb_degree_of_freedom * nb_nodes_per_element);

      fem_boundary.computeNtb(dual_before_integ, dual_by_shapes, type,
                              ghost_type, element_ids);

      Array<Real> dual_by_shapes_integ(nb_elements, nb_degree_of_freedom *
                                                        nb_nodes_per_element);
      fem_boundary.integrate(dual_by_shapes, dual_by_shapes_integ,
                             nb_degree_of_freedom * nb_nodes_per_element, type,
                             ghost_type, element_ids);

      // assemble the result into force vector
      model.getDOFManager().assembleElementalArrayLocalArray(
          dual_by_shapes_integ, dual, type, ghost_type, 1., element_ids);
    }
  }
};

/* -------------------------------------------------------------------------- */
template <typename ModelType>
template <typename FunctorType>
inline void BoundaryCondition<ModelType>::applyBC(FunctorType && func) {
  auto bit = model->getMesh().getGroupManager().element_group_begin();
  auto bend = model->getMesh().getGroupManager().element_group_end();
  for (; bit != bend; ++bit) {
    applyBC(std::forward<decltype(func)>(func), *bit);
  }
}

/* -------------------------------------------------------------------------- */
template <typename ModelType>
template <typename FunctorType>
inline void
BoundaryCondition<ModelType>::applyBC(FunctorType && func,
                                      const std::string & group_name) {
  try {
    const ElementGroup & element_group =
        model->getMesh().getElementGroup(group_name);
    applyBC(std::forward<decltype(func)>(func), element_group);
  } catch (akantu::debug::Exception & e) {
    AKANTU_EXCEPTION("Error applying a boundary condition onto \""
                     << group_name << "\"! [" << e.what() << "]");
  }
}

/* -------------------------------------------------------------------------- */
template <typename ModelType>
template <typename FunctorType>
inline void
BoundaryCondition<ModelType>::applyBC(FunctorType && func,
                                      const ElementGroup & element_group) {
#if !defined(AKANTU_NDEBUG)
  if (element_group.getDimension() != model->getSpatialDimension() - 1) {
    AKANTU_DEBUG_INFO("The group "
                      << element_group.getName()
                      << " does not contain only boundaries elements");
  }
#endif

  TemplateFunctionWrapper<FunctorType>::applyBC(
      std::forward<decltype(func)>(func), element_group, *this);
}

#endif /* AKANTU_BOUNDARY_CONDITION_TMPL_HH_ */

} // namespace akantu
