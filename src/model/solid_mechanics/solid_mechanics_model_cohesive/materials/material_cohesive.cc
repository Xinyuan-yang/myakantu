/**
 * Copyright (©) 2012-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "material_cohesive.hh"
#include "aka_random_generator.hh"
#include "dof_synchronizer.hh"
#include "fe_engine_template.hh"
#include "integrator_gauss.hh"
#include "shape_cohesive.hh"
#include "sparse_matrix.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
MaterialCohesive::MaterialCohesive(SolidMechanicsModel & model, const ID & id)
    : Material(model, id, "SolidMechanicsFEEngine"),
      facet_filter(
          std::make_shared<ElementTypeMapArray<Idx>>("facet_filter", id)),
      fem_cohesive(
          model.getFEEngineClass<MyFEEngineCohesiveType>("CohesiveFEEngine")),
      reversible_energy(registerInternal<Real, CohesiveInternalField>(
          "reversible_energy", 1)),
      total_energy(
          registerInternal<Real, CohesiveInternalField>("total_energy", 1)),
      opening(registerInternal<Real, CohesiveInternalField>("opening",
                                                            spatial_dimension)),
      tractions(registerInternal<Real, CohesiveInternalField>(
          "tractions", spatial_dimension)),
      contact_tractions(registerInternal<Real, CohesiveInternalField>(
          "contact_tractions", spatial_dimension)),
      contact_opening(registerInternal<Real, CohesiveInternalField>(
          "contact_opening", spatial_dimension)),
      delta_max(registerInternal<Real, CohesiveInternalField>("delta max", 1)),
      use_previous_delta_max(false), use_previous_opening(false),
      damage(registerInternal<Real, CohesiveInternalField>("damage", 1)),
      sigma_c(registerInternal<Real, FacetRandomInternalField>("sigma_c", 1)),
      model(&aka::as_type<SolidMechanicsModelCohesive>(model)),
      normals(registerInternal<Real, CohesiveInternalField>(
          "normal", spatial_dimension)) {
  AKANTU_DEBUG_IN();

  this->registerParam("sigma_c", sigma_c, _pat_parsable | _pat_readable,
                      "Critical stress");
  this->registerParam("delta_c", delta_c, Real(0.),
                      _pat_parsable | _pat_readable, "Critical displacement");

  this->getElementFilter().initialize(this->model->getMesh(),
                                      _spatial_dimension = spatial_dimension,
                                      _element_kind = _ek_cohesive);

  if (this->model->getIsExtrinsic()) {
    this->getFacetFilter().initialize(this->model->getMeshFacets(),
                                      _spatial_dimension =
                                          spatial_dimension - 1,
                                      _element_kind = _ek_regular);
  }

  if (this->use_previous_delta_max) {
    this->delta_max.initializeHistory();
  }

  this->opening.initializeHistory();
  this->tractions.initializeHistory();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void MaterialCohesive::assembleInternalForces(GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  auto & internal_force = const_cast<Array<Real> &>(model->getInternalForce());

  for (auto type : getElementFilter().elementTypes(spatial_dimension,
                                                   ghost_type, _ek_cohesive)) {
    auto & elem_filter = getElementFilter(type, ghost_type);
    auto nb_element = elem_filter.size();
    if (nb_element == 0) {
      continue;
    }

    const auto & shapes = fem_cohesive.getShapes(type, ghost_type);
    auto & traction = tractions(type, ghost_type);
    auto & contact_traction = contact_tractions(type, ghost_type);

    auto size_of_shapes = shapes.getNbComponent();
    auto nb_nodes_per_element = Mesh::getNbNodesPerElement(type);
    auto nb_quadrature_points =
        fem_cohesive.getNbIntegrationPoints(type, ghost_type);

    /// compute @f$t_i N_a@f$

    auto traction_cpy = std::make_shared<Array<Real>>(
        nb_element * nb_quadrature_points, spatial_dimension * size_of_shapes);

    auto traction_it = traction.begin(spatial_dimension, 1);
    auto contact_traction_it = contact_traction.begin(spatial_dimension, 1);
    auto shapes_filtered_begin = shapes.begin(1, size_of_shapes);
    auto traction_cpy_it =
        traction_cpy->begin(spatial_dimension, size_of_shapes);

    Matrix<Real> traction_tmp(spatial_dimension, 1);

    for (Int el = 0; el < nb_element; ++el) {
      auto current_quad = elem_filter(el) * nb_quadrature_points;

      for (Int q = 0; q < nb_quadrature_points; ++q, ++traction_it,
               ++contact_traction_it, ++current_quad, ++traction_cpy_it) {

        auto && shapes_filtered = shapes_filtered_begin[current_quad];

        *traction_cpy_it =
            (*traction_it + *contact_traction_it) * shapes_filtered;
      }
    }

    /**
     * compute @f$\int t \cdot N\, dS@f$ by  @f$ \sum_q \mathbf{N}^t
     * \mathbf{t}_q \overline w_q J_q@f$
     */
    auto partial_int_t_N = std::make_shared<Array<Real>>(
        nb_element, spatial_dimension * size_of_shapes, "int_t_N");

    fem_cohesive.integrate(*traction_cpy, *partial_int_t_N,
                           spatial_dimension * size_of_shapes, type, ghost_type,
                           elem_filter);

    auto int_t_N = std::make_shared<Array<Real>>(
        nb_element, 2 * spatial_dimension * size_of_shapes, "int_t_N");

    auto * int_t_N_val = int_t_N->data();
    auto * partial_int_t_N_val = partial_int_t_N->data();
    for (Int el = 0; el < nb_element; ++el) {
      std::copy_n(partial_int_t_N_val, size_of_shapes * spatial_dimension,
                  int_t_N_val);
      std::copy_n(partial_int_t_N_val, size_of_shapes * spatial_dimension,
                  int_t_N_val + size_of_shapes * spatial_dimension);

      for (Int n = 0; n < size_of_shapes * spatial_dimension; ++n) {
        int_t_N_val[n] *= -1.;
      }

      int_t_N_val += nb_nodes_per_element * spatial_dimension;
      partial_int_t_N_val += size_of_shapes * spatial_dimension;
    }

    /// assemble
    model->getDOFManager().assembleElementalArrayLocalArray(
        *int_t_N, internal_force, type, ghost_type, 1, elem_filter);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void MaterialCohesive::assembleStiffnessMatrix(GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  for (auto type : getElementFilter().elementTypes(spatial_dimension,
                                                   ghost_type, _ek_cohesive)) {
    auto nb_quadrature_points =
        fem_cohesive.getNbIntegrationPoints(type, ghost_type);
    auto nb_nodes_per_element = Mesh::getNbNodesPerElement(type);

    auto & elem_filter = getElementFilter(type, ghost_type);
    auto nb_element = elem_filter.size();

    if (nb_element == 0U) {
      continue;
    }

    const auto & shapes = fem_cohesive.getShapes(type, ghost_type);
    auto size_of_shapes = shapes.getNbComponent();

    auto shapes_filtered = std::make_shared<Array<Real>>(
        nb_element * nb_quadrature_points, size_of_shapes, "filtered shapes");

    for (auto && data :
         zip(filter(elem_filter,
                    make_view(shapes, size_of_shapes, nb_quadrature_points)),
             make_view(*shapes_filtered, size_of_shapes,
                       nb_quadrature_points))) {
      std::get<1>(data) = std::get<0>(data);
    }

    Matrix<Real> A(spatial_dimension * size_of_shapes,
                   spatial_dimension * nb_nodes_per_element);
    A.zero();
    for (Int i = 0; i < spatial_dimension * size_of_shapes; ++i) {
      A(i, i) = 1;
      A(i, i + spatial_dimension * size_of_shapes) = -1;
    }

    /// get the tangent matrix @f$\frac{\partial{(t/\delta)}}{\partial{\delta}}
    /// @f$
    auto tangent_stiffness_matrix = std::make_unique<Array<Real>>(
        nb_element * nb_quadrature_points,
        spatial_dimension * spatial_dimension, "tangent_stiffness_matrix");

    computeNormal(model->getCurrentPosition(), normals(type, ghost_type), type,
                  ghost_type);

    /// compute openings @f$\mathbf{\delta}@f$
    // computeOpening(model->getDisplacement(), opening(type, ghost_type), type,
    // ghost_type);

    tangent_stiffness_matrix->zero();

    computeTangentTraction(type, *tangent_stiffness_matrix, ghost_type);

    UInt size_at_nt_d_n_a = spatial_dimension * nb_nodes_per_element *
                            spatial_dimension * nb_nodes_per_element;
    auto at_nt_d_n_a = std::make_unique<Array<Real>>(
        nb_element * nb_quadrature_points, size_at_nt_d_n_a, "A^t*N^t*D*N*A");

    Matrix<Real> N(spatial_dimension, spatial_dimension * size_of_shapes);

    for (auto && [At_Nt_D_N_A, D, shapes] :
         zip(make_view(*at_nt_d_n_a, spatial_dimension * nb_nodes_per_element,
                       spatial_dimension * nb_nodes_per_element),
             make_view(*tangent_stiffness_matrix, spatial_dimension,
                       spatial_dimension),
             make_view(*shapes_filtered, size_of_shapes))) {
      N.zero();
      /**
       * store  the   shapes  in  voigt   notations  matrix  @f$\mathbf{N}  =
       * \begin{array}{cccccc} N_0(\xi) & 0 & N_1(\xi)  &0 & N_2(\xi) & 0 \\
       * 0 & * N_0(\xi)& 0 &N_1(\xi)& 0 & N_2(\xi) \end{array} @f$
       **/
      for (Int i = 0; i < spatial_dimension; ++i) {
        for (Int n = 0; n < size_of_shapes; ++n) {
          N(i, i + spatial_dimension * n) = shapes(n);
        }
      }

      /**
       * compute stiffness matrix  @f$   \mathbf{K}    =    \delta
       *\mathbf{U}^T \int_{\Gamma_c}    {\mathbf{P}^t
       *\frac{\partial{\mathbf{t}}}
       *{\partial{\delta}}
       * \mathbf{P} d\Gamma \Delta \mathbf{U}}  @f$
       **/
      Matrix<Real> NA = N * A;
      At_Nt_D_N_A = (D * NA).transpose() * NA;
    }

    auto K_e =
        std::make_unique<Array<Real>>(nb_element, size_at_nt_d_n_a, "K_e");

    fem_cohesive.integrate(*at_nt_d_n_a, *K_e, size_at_nt_d_n_a, type,
                           ghost_type, elem_filter);

    model->getDOFManager().assembleElementalMatricesToMatrix(
        "K", "displacement", *K_e, type, ghost_type, _unsymmetric, elem_filter);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- *
 * Compute traction from displacements
 *
 * @param[in] ghost_type compute the residual for _ghost or _not_ghost element
 */
void MaterialCohesive::computeTraction(GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  for (const auto & type : getElementFilter().elementTypes(
           spatial_dimension, ghost_type, _ek_cohesive)) {
    auto & elem_filter = getElementFilter(type, ghost_type);
    auto nb_element = elem_filter.size();
    if (nb_element == 0) {
      continue;
    }

    /// compute normals @f$\mathbf{n}@f$
    computeNormal(model->getCurrentPosition(), normals(type, ghost_type), type,
                  ghost_type);

    /// compute openings @f$\mathbf{\delta}@f$
    computeOpening(model->getDisplacement(), opening(type, ghost_type), type,
                   ghost_type);

    /// compute traction @f$\mathbf{t}@f$
    computeTraction(type, ghost_type);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void MaterialCohesive::computeNormal(const Array<Real> & position,
                                     Array<Real> & normal, ElementType type,
                                     GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  auto & fem_cohesive =
      this->model->getFEEngineClass<MyFEEngineCohesiveType>("CohesiveFEEngine");

  normal.zero();

  tuple_dispatch<ElementTypes_t<_ek_cohesive>>(
      [&](auto && enum_type) {
        constexpr ElementType type = aka::decay_v<decltype(enum_type)>;
        fem_cohesive.getShapeFunctions()
            .computeNormalsOnIntegrationPoints<type,
                                               CohesiveReduceFunctionMean>(
                position, normal, ghost_type,
                getElementFilter(type, ghost_type));
      },
      type);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void MaterialCohesive::computeOpening(const Array<Real> & displacement,
                                      Array<Real> & opening, ElementType type,
                                      GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  auto & fem_cohesive =
      this->model->getFEEngineClass<MyFEEngineCohesiveType>("CohesiveFEEngine");

  tuple_dispatch<ElementTypes_t<_ek_cohesive>>(
      [&](auto && enum_type) {
        constexpr ElementType type = aka::decay_v<decltype(enum_type)>;
        fem_cohesive.getShapeFunctions()
            .interpolateOnIntegrationPoints<type,
                                            CohesiveReduceFunctionOpening>(
                displacement, opening, spatial_dimension, ghost_type,
                getElementFilter(type, ghost_type));
      },
      type);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void MaterialCohesive::updateEnergies(ElementType type) {
  AKANTU_DEBUG_IN();

  if (Mesh::getKind(type) != _ek_cohesive) {
    return;
  }

  Vector<Real> b(spatial_dimension);
  Vector<Real> h(spatial_dimension);
  auto erev = reversible_energy(type).begin();
  auto etot = total_energy(type).begin();
  auto traction_it = tractions(type).begin(spatial_dimension);
  auto traction_old_it = tractions.previous(type).begin(spatial_dimension);
  auto opening_it = opening(type).begin(spatial_dimension);
  auto opening_old_it = opening.previous(type).begin(spatial_dimension);

  auto traction_end = tractions(type).end(spatial_dimension);

  /// loop on each quadrature point
  for (; traction_it != traction_end; ++traction_it, ++traction_old_it,
                                      ++opening_it, ++opening_old_it, ++erev,
                                      ++etot) {
    /// trapezoidal integration
    b = *opening_it;
    b -= *opening_old_it;

    h = *traction_old_it;
    h += *traction_it;

    *etot += .5 * b.dot(h);
    *erev = .5 * traction_it->dot(*opening_it);
  }

  /// update old values

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
Real MaterialCohesive::getReversibleEnergy() {
  AKANTU_DEBUG_IN();
  Real erev = 0.;

  /// integrate reversible energy for each type of elements
  for (const auto & type : getElementFilter().elementTypes(
           spatial_dimension, _not_ghost, _ek_cohesive)) {
    erev +=
        fem_cohesive.integrate(reversible_energy(type, _not_ghost), type,
                               _not_ghost, getElementFilter(type, _not_ghost));
  }

  AKANTU_DEBUG_OUT();
  return erev;
}

/* -------------------------------------------------------------------------- */
Real MaterialCohesive::getDissipatedEnergy() {
  AKANTU_DEBUG_IN();
  Real edis = 0.;

  /// integrate dissipated energy for each type of elements
  for (const auto & type : getElementFilter().elementTypes(
           spatial_dimension, _not_ghost, _ek_cohesive)) {
    Array<Real> dissipated_energy(total_energy(type, _not_ghost));
    dissipated_energy -= reversible_energy(type, _not_ghost);
    edis += fem_cohesive.integrate(dissipated_energy, type, _not_ghost,
                                   getElementFilter(type, _not_ghost));
  }

  AKANTU_DEBUG_OUT();
  return edis;
}

/* -------------------------------------------------------------------------- */
Real MaterialCohesive::getContactEnergy() {
  AKANTU_DEBUG_IN();
  Real econ = 0.;

  /// integrate contact energy for each type of elements
  for (const auto & type : getElementFilter().elementTypes(
           spatial_dimension, _not_ghost, _ek_cohesive)) {

    auto & el_filter = getElementFilter(type, _not_ghost);
    auto nb_quad_per_el = fem_cohesive.getNbIntegrationPoints(type, _not_ghost);
    auto nb_quad_points = el_filter.size() * nb_quad_per_el;
    Array<Real> contact_energy(nb_quad_points);

    auto contact_traction_it =
        contact_tractions(type, _not_ghost).begin(spatial_dimension);
    auto contact_opening_it =
        contact_opening(type, _not_ghost).begin(spatial_dimension);

    /// loop on each quadrature point
    for (Int q = 0; q < nb_quad_points;
         ++contact_traction_it, ++contact_opening_it, ++q) {

      contact_energy(q) = .5 * contact_traction_it->dot(*contact_opening_it);
    }

    econ += fem_cohesive.integrate(contact_energy, type, _not_ghost, el_filter);
  }

  AKANTU_DEBUG_OUT();
  return econ;
}

/* -------------------------------------------------------------------------- */
Real MaterialCohesive::getEnergy(const std::string & type) {
  if (type == "reversible") {
    return getReversibleEnergy();
  }
  if (type == "dissipated") {
    return getDissipatedEnergy();
  }
  if (type == "cohesive contact") {
    return getContactEnergy();
  }

  return 0.;
}

/* -------------------------------------------------------------------------- */
} // namespace akantu
