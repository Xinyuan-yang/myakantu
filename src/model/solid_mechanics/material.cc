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
#include "material.hh"
#include "mesh_iterators.hh"
#include "solid_mechanics_model.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
Material::Material(SolidMechanicsModel & model, const ID & id,
                   const ID & fe_engine_id)
    : Parent(model, id, model.getSpatialDimension(), _ek_regular, fe_engine_id),
      stress(registerInternal("stress", spatial_dimension * spatial_dimension,
                              fe_engine_id)),
      eigengradu(registerInternal(
          "eigen_grad_u", spatial_dimension * spatial_dimension, fe_engine_id)),
      gradu(registerInternal("grad_u", spatial_dimension * spatial_dimension,
                             fe_engine_id)),
      potential_energy(registerInternal("potential_energy", 1, fe_engine_id)),
      interpolation_inverse_coordinates("interpolation inverse coordinates",
                                        id),
      interpolation_points_matrices("interpolation points matrices", id),
      eigen_grad_u(model.getSpatialDimension(), model.getSpatialDimension()) {
  eigen_grad_u.setZero();
  registerParam("rho", rho, Real(0.), _pat_parsable | _pat_modifiable,
                "Density");
  registerParam("finite_deformation", finite_deformation, false,
                _pat_parsable | _pat_readable, "Is finite deformation");
  registerParam("inelastic_deformation", inelastic_deformation, false,
                _pat_internal, "Is inelastic deformation");
  registerParam("eigen_grad_u", eigen_grad_u, _pat_parsable, "EigenGradU");

  this->getModel().registerEventHandler(*this);
}

/* -------------------------------------------------------------------------- */
void Material::initMaterial() {
  AKANTU_DEBUG_IN();

  if (finite_deformation) {
    this->registerInternal("piola_kirchhoff_2",
                           spatial_dimension * spatial_dimension);
    this->piola_kirchhoff_2 = this->getSharedPtrInternal("piola_kirchhoff_2");
    this->piola_kirchhoff_2->initializeHistory();
    this->registerInternal("green_strain",
                           spatial_dimension * spatial_dimension);
    this->green_strain = this->getSharedPtrInternal("green_strain");
  }

  this->stress.initializeHistory();
  this->gradu.initializeHistory();

  Parent::initConstitutiveLaw();
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void Material::updateInternalParameters() {
  auto dim = getModel().getSpatialDimension();
  for (const auto & type :
       getElementFilter().elementTypes(_element_kind = _ek_regular)) {
    for (auto eigen_gradu : make_view(eigengradu(type), dim, dim)) {
      eigen_gradu = eigen_grad_u;
    }
  }

  Parent::updateInternalParameters();
}

/* -------------------------------------------------------------------------- */
/**
 * Compute the internal forces by assembling @f$\int_{e} \sigma_e \frac{\partial
 * \varphi}{\partial X} dX @f$
 *
 * @param[in] ghost_type compute the internal forces for _ghost or _not_ghost
 * element
 */
void Material::assembleInternalForces(GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  Int spatial_dimension = getModel().getSpatialDimension();

  tuple_dispatch<AllSpatialDimensions>(
      [&](auto && _) {
        constexpr auto dim = aka::decay_v<decltype(_)>;

        for (auto && type :
             getElementFilter().elementTypes(spatial_dimension, ghost_type)) {
          if (not finite_deformation) {
            this->assembleInternalForces<dim>(type, ghost_type);
          } else {
            this->assembleInternalForcesFiniteDeformation<dim>(type,
                                                               ghost_type);
          }
        }
      },
      spatial_dimension);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
/**
 * Compute  the  stress from the gradu
 *
 * @param[in] ghost_type compute the residual for _ghost or _not_ghost element
 */
void Material::computeAllStresses(GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  auto spatial_dimension = getModel().getSpatialDimension();
  auto & fem = getFEEngine();
  for (const auto & type :
       getElementFilter().elementTypes(spatial_dimension, ghost_type)) {
    auto & elem_filter = getElementFilter(type, ghost_type);

    if (elem_filter.empty()) {
      continue;
    }
    auto & gradu_vect = gradu(type, ghost_type);

    /// compute @f$\nabla u@f$
    fem.gradientOnIntegrationPoints(getModel().getDisplacement(), gradu_vect,
                                    spatial_dimension, type, ghost_type,
                                    elem_filter);

    gradu_vect -= eigengradu(type, ghost_type);

    /// compute @f$\mathbf{\sigma}_q@f$ from @f$\nabla u@f$
    computeStress(type, ghost_type);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void Material::computeAllCauchyStresses(GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_ASSERT(finite_deformation, "The Cauchy stress can only be "
                                          "computed if you are working in "
                                          "finite deformation.");

  for (auto type :
       getElementFilter().elementTypes(spatial_dimension, ghost_type)) {
    tuple_dispatch<AllSpatialDimensions>(
        [&](auto && _) {
          constexpr auto dim = aka::decay_v<decltype(_)>;
          this->StoCauchy<dim>(type, ghost_type);
        },
        spatial_dimension);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <Int dim>
void Material::StoCauchy(ElementType el_type, GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  for (auto && [grad_u, piola, sigma] :
       zip(make_view<dim, dim>(this->gradu(el_type, ghost_type)),
           make_view<dim, dim>((*this->piola_kirchhoff_2)(el_type, ghost_type)),
           make_view<dim, dim>(this->stress(el_type, ghost_type)))) {
    Matrix<Real, dim, dim> F = gradUToF<dim>(grad_u);
    this->StoCauchy<dim>(F, piola, sigma);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void Material::setToSteadyState(GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  const auto & displacement = getModel().getDisplacement();
  auto & fem = getFEEngine();
  // resizeInternalArray(gradu);

  auto spatial_dimension = getModel().getSpatialDimension();

  for (auto type :
       getElementFilter().elementTypes(spatial_dimension, ghost_type)) {
    auto & elem_filter = getElementFilter(type, ghost_type);
    auto & gradu_vect = gradu(type, ghost_type);

    /// compute @f$\nabla u@f$
    fem.gradientOnIntegrationPoints(displacement, gradu_vect, spatial_dimension,
                                    type, ghost_type, elem_filter);

    setToSteadyState(type, ghost_type);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
/**
 * Compute  the stiffness  matrix by  assembling @f$\int_{\omega}  B^t  \times D
 * \times B d\omega @f$
 *
 * @param[in] ghost_type compute the residual for _ghost or _not_ghost element
 */
void Material::assembleStiffnessMatrix(GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  auto spatial_dimension = getModel().getSpatialDimension();

  tuple_dispatch<AllSpatialDimensions>(
      [&](auto && _) {
        constexpr auto dim = aka::decay_v<decltype(_)>;

        for (auto type :
             getElementFilter().elementTypes(spatial_dimension, ghost_type)) {
          if (finite_deformation) {
            this->assembleStiffnessMatrixFiniteDeformation<dim>(type,
                                                                ghost_type);
          } else {
            this->assembleStiffnessMatrix<dim>(type, ghost_type);
          }
        }
      },
      spatial_dimension);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <Int dim>
void Material::assembleStiffnessMatrix(ElementType type, GhostType ghost_type) {
  tuple_dispatch<AllElementTypes>(
      [&](auto && enum_type) {
        constexpr auto type = aka::decay_v<decltype(enum_type)>;
        this->assembleStiffnessMatrix<dim, type>(ghost_type);
      },
      type);
}

/* -------------------------------------------------------------------------- */
template <Int dim, ElementType type>
void Material::assembleStiffnessMatrix(GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  const auto & elem_filter = getElementFilter(type, ghost_type);
  if (elem_filter.empty()) {
    AKANTU_DEBUG_OUT();
    return;
  }

  auto & fem = getFEEngine();
  auto & gradu_vect = gradu(type, ghost_type);

  auto nb_element = elem_filter.size();
  auto nb_quadrature_points = fem.getNbIntegrationPoints(type, ghost_type);
  constexpr auto nb_nodes_per_element = Mesh::getNbNodesPerElement(type);
  constexpr auto tangent_size = getTangentStiffnessVoigtSize(dim);

  gradu_vect.resize(nb_quadrature_points * nb_element);

  fem.gradientOnIntegrationPoints(getModel().getDisplacement(), gradu_vect, dim,
                                  type, ghost_type, elem_filter);

  auto tangent_stiffness_matrix = std::make_unique<Array<Real>>(
      nb_element * nb_quadrature_points, tangent_size * tangent_size,
      "tangent_stiffness_matrix");

  tangent_stiffness_matrix->zero();

  computeTangentModuli(type, *tangent_stiffness_matrix, ghost_type);

  /// compute @f$\mathbf{B}^t * \mathbf{D} * \mathbf{B}@f$
  auto bt_d_b_size = dim * nb_nodes_per_element;

  auto bt_d_b = std::make_unique<Array<Real>>(
      nb_element * nb_quadrature_points, bt_d_b_size * bt_d_b_size, "B^t*D*B");

  fem.computeBtDB(*tangent_stiffness_matrix, *bt_d_b, 4, type, ghost_type,
                  elem_filter);

  /// compute @f$ k_e = \int_e \mathbf{B}^t * \mathbf{D} * \mathbf{B}@f$
  auto K_e = std::make_unique<Array<Real>>(nb_element,
                                           bt_d_b_size * bt_d_b_size, "K_e");

  fem.integrate(*bt_d_b, *K_e, bt_d_b_size * bt_d_b_size, type, ghost_type,
                elem_filter);

  getModel().getDOFManager().assembleElementalMatricesToMatrix(
      "K", "displacement", *K_e, type, ghost_type, _symmetric, elem_filter);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <Int dim>
void Material::assembleStiffnessMatrixFiniteDeformation(ElementType type,
                                                        GhostType ghost_type) {
  tuple_dispatch<AllElementTypes>(
      [&](auto && enum_type) {
        constexpr auto type = aka::decay_v<decltype(enum_type)>;
        this->assembleStiffnessMatrixNL<dim, type>(ghost_type);
        this->assembleStiffnessMatrixL2<dim, type>(ghost_type);
      },
      type);
}

/* -------------------------------------------------------------------------- */
template <Int dim, ElementType type>
void Material::assembleStiffnessMatrixNL(GhostType ghost_type) {
  AKANTU_DEBUG_IN();
  auto & fem = getFEEngine();

  const auto & shapes_derivatives = fem.getShapesDerivatives(type, ghost_type);
  const auto & elem_filter = getElementFilter(type, ghost_type);

  const auto nb_element = elem_filter.size();
  constexpr auto nb_nodes_per_element = Mesh::getNbNodesPerElement(type);
  const auto nb_quadrature_points =
      fem.getNbIntegrationPoints(type, ghost_type);

  auto shapes_derivatives_filtered = std::make_unique<Array<Real>>(
      nb_element * nb_quadrature_points, dim * nb_nodes_per_element,
      "shapes derivatives filtered");

  FEEngine::filterElementalData(fem.getMesh(), shapes_derivatives,
                                *shapes_derivatives_filtered, type, ghost_type,
                                elem_filter);

  /// compute @f$\mathbf{B}^t * \mathbf{D} * \mathbf{B}@f$
  constexpr auto bt_s_b_size = dim * nb_nodes_per_element;
  constexpr auto piola_matrix_size = getCauchyStressMatrixSize(dim);

  auto bt_s_b = std::make_unique<Array<Real>>(
      nb_element * nb_quadrature_points, bt_s_b_size * bt_s_b_size, "B^t*D*B");

  Matrix<Real, piola_matrix_size, bt_s_b_size> B;
  Matrix<Real, piola_matrix_size, piola_matrix_size> S;

  for (auto && [Bt_S_B, Piola_kirchhoff_matrix, shapes_derivatives] :
       zip(make_view<bt_s_b_size, bt_s_b_size>(*bt_s_b),
           make_view<dim, dim>((*piola_kirchhoff_2)(type, ghost_type)),
           make_view<dim, nb_nodes_per_element>(
               *shapes_derivatives_filtered))) {
    setCauchyStressMatrix<dim>(Piola_kirchhoff_matrix, S);
    VoigtHelper<dim>::transferBMatrixToBNL(shapes_derivatives, B,
                                           nb_nodes_per_element);
    Bt_S_B = B.transpose() * S * B;
  }

  /// compute @f$ k_e = \int_e \mathbf{B}^t * \mathbf{D} * \mathbf{B}@f$
  auto K_e = std::make_unique<Array<Real>>(nb_element,
                                           bt_s_b_size * bt_s_b_size, "K_e");

  fem.integrate(*bt_s_b, *K_e, bt_s_b_size * bt_s_b_size, type, ghost_type,
                elem_filter);

  getModel().getDOFManager().assembleElementalMatricesToMatrix(
      "K", "displacement", *K_e, type, ghost_type, _symmetric, elem_filter);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <Int dim, ElementType type>
void Material::assembleStiffnessMatrixL2(GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  auto & fem = getFEEngine();
  const auto & shapes_derivatives = fem.getShapesDerivatives(type, ghost_type);

  auto & elem_filter = getElementFilter(type, ghost_type);
  auto & gradu_vect = gradu(type, ghost_type);

  auto nb_element = elem_filter.size();
  constexpr auto nb_nodes_per_element = Mesh::getNbNodesPerElement(type);
  auto nb_quadrature_points = fem.getNbIntegrationPoints(type, ghost_type);

  gradu_vect.resize(nb_quadrature_points * nb_element);

  fem.gradientOnIntegrationPoints(getModel().getDisplacement(), gradu_vect, dim,
                                  type, ghost_type, elem_filter);

  constexpr auto tangent_size = getTangentStiffnessVoigtSize(dim);

  auto tangent_stiffness_matrix = std::make_unique<Array<Real>>(
      nb_element * nb_quadrature_points, tangent_size * tangent_size,
      "tangent_stiffness_matrix");

  tangent_stiffness_matrix->zero();

  computeTangentModuli(type, *tangent_stiffness_matrix, ghost_type);

  auto shapes_derivatives_filtered = std::make_unique<Array<Real>>(
      nb_element * nb_quadrature_points, dim * nb_nodes_per_element,
      "shapes derivatives filtered");
  FEEngine::filterElementalData(fem.getMesh(), shapes_derivatives,
                                *shapes_derivatives_filtered, type, ghost_type,
                                elem_filter);

  /// compute @f$\mathbf{B}^t * \mathbf{D} * \mathbf{B}@f$
  constexpr auto bt_d_b_size = dim * nb_nodes_per_element;
  auto bt_d_b = std::make_unique<Array<Real>>(
      nb_element * nb_quadrature_points, bt_d_b_size * bt_d_b_size, "B^t*D*B");

  Matrix<Real, tangent_size, bt_d_b_size> B;
  Matrix<Real, tangent_size, bt_d_b_size> B2;

  for (auto && [Bt_D_B, grad_u, D, shapes_derivative] :
       zip(make_view<bt_d_b_size, bt_d_b_size>(*bt_d_b),
           make_view<dim, dim>(gradu_vect),
           make_view<tangent_size, tangent_size>(*tangent_stiffness_matrix),
           make_view<dim, nb_nodes_per_element>(
               *shapes_derivatives_filtered))) {
    VoigtHelper<dim>::transferBMatrixToSymVoigtBMatrix(shapes_derivative, B,
                                                       nb_nodes_per_element);
    VoigtHelper<dim>::transferBMatrixToBL2(shapes_derivative, grad_u, B2,
                                           nb_nodes_per_element);
    B += B2;
    Bt_D_B = B.transpose() * D * B;
  }

  /// compute @f$ k_e = \int_e \mathbf{B}^t * \mathbf{D} * \mathbf{B}@f$
  auto K_e = std::make_unique<Array<Real>>(nb_element,
                                           bt_d_b_size * bt_d_b_size, "K_e");

  fem.integrate(*bt_d_b, *K_e, bt_d_b_size * bt_d_b_size, type, ghost_type,
                elem_filter);

  getModel().getDOFManager().assembleElementalMatricesToMatrix(
      "K", "displacement", *K_e, type, ghost_type, _symmetric, elem_filter);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <Int dim>
void Material::assembleInternalForces(ElementType type, GhostType ghost_type) {
  tuple_dispatch<AllElementTypes>(
      [&](auto && enum_type) {
        constexpr auto type = aka::decay_v<decltype(enum_type)>;
        this->assembleInternalForces<dim, type>(ghost_type);
      },
      type);
}

/* -------------------------------------------------------------------------- */
template <Int dim, ElementType type>
void Material::assembleInternalForces(GhostType ghost_type) {
  auto & internal_force = getModel().getInternalForce();

  // Mesh & mesh = fem.getMesh();
  auto && elem_filter = getElementFilter(type, ghost_type);
  auto nb_element = elem_filter.size();

  if (nb_element == 0) {
    return;
  }

  auto && fem = this->getFEEngine();
  const auto & shapes_derivatives = fem.getShapesDerivatives(type, ghost_type);

  auto size_of_shapes_derivatives = shapes_derivatives.getNbComponent();
  auto nb_quadrature_points = fem.getNbIntegrationPoints(type, ghost_type);
  auto nb_nodes_per_element = Mesh::getNbNodesPerElement(type);

  /// compute @f$\sigma \frac{\partial \varphi}{\partial X}@f$ by
  /// @f$\mathbf{B}^t \mathbf{\sigma}_q@f$
  auto sigma_dphi_dx = std::make_unique<Array<Real>>(
      nb_element * nb_quadrature_points, size_of_shapes_derivatives,
      "sigma_x_dphi_/_dX");

  fem.computeBtD(stress(type, ghost_type), *sigma_dphi_dx, type, ghost_type,
                 elem_filter);

  /**
   * compute @f$\int \sigma  * \frac{\partial \varphi}{\partial X}dX@f$ by
   * @f$ \sum_q \mathbf{B}^t
   * \mathbf{\sigma}_q \overline w_q J_q@f$
   */
  auto int_sigma_dphi_dx = std::make_unique<Array<Real>>(
      nb_element, nb_nodes_per_element * spatial_dimension,
      "int_sigma_x_dphi_/_dX");

  fem.integrate(*sigma_dphi_dx, *int_sigma_dphi_dx, size_of_shapes_derivatives,
                type, ghost_type, elem_filter);

  /// assemble
  getModel().getDOFManager().assembleElementalArrayLocalArray(
      *int_sigma_dphi_dx, internal_force, type, ghost_type, -1, elem_filter);
}
/* -------------------------------------------------------------------------- */
template <Int dim>
void Material::assembleInternalForcesFiniteDeformation(ElementType type,
                                                       GhostType ghost_type) {
  tuple_dispatch<AllElementTypes>(
      [&](auto && enum_type) {
        constexpr auto type = aka::decay_v<decltype(enum_type)>;
        this->assembleInternalForcesFiniteDeformation<dim, type>(ghost_type);
      },
      type);
}

/* -------------------------------------------------------------------------- */
template <Int dim, ElementType type>
void Material::assembleInternalForcesFiniteDeformation(GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  auto & fem = getFEEngine();
  auto & internal_force = getModel().getInternalForce();
  auto & mesh = fem.getMesh();
  const auto & shapes_derivatives = fem.getShapesDerivatives(type, ghost_type);

  auto & elem_filter = getElementFilter(type, ghost_type);
  if (elem_filter.empty()) {
    return;
  }

  auto size_of_shapes_derivatives = shapes_derivatives.getNbComponent();
  auto nb_element = elem_filter.size();

  constexpr auto nb_nodes_per_element = Mesh::getNbNodesPerElement(type);
  constexpr auto stress_size = getTangentStiffnessVoigtSize(dim);
  constexpr auto bt_s_size = dim * nb_nodes_per_element;

  auto nb_quadrature_points = fem.getNbIntegrationPoints(type, ghost_type);

  auto shapesd_filtered = std::make_unique<Array<Real>>(
      nb_element, size_of_shapes_derivatives, "filtered shapesd");

  FEEngine::filterElementalData(mesh, shapes_derivatives, *shapesd_filtered,
                                type, ghost_type, elem_filter);

  // Set stress vectors
  auto bt_s = std::make_unique<Array<Real>>(nb_element * nb_quadrature_points,
                                            bt_s_size, "B^t*S");

  Matrix<Real, stress_size, bt_s_size> B_tensor;
  Matrix<Real, stress_size, bt_s_size> B2_tensor;

  for (auto && [grad_u, S, r, shapes_derivative] :
       zip(make_view<dim, dim>(this->gradu(type, ghost_type)),
           make_view<dim, dim>((*this->piola_kirchhoff_2)(type, ghost_type)),
           make_view<bt_s_size>(*bt_s),
           make_view<dim, nb_nodes_per_element>(*shapesd_filtered))) {
    VoigtHelper<dim>::transferBMatrixToSymVoigtBMatrix(
        shapes_derivative, B_tensor, nb_nodes_per_element);

    VoigtHelper<dim>::transferBMatrixToBL2(shapes_derivative, grad_u, B2_tensor,
                                           nb_nodes_per_element);

    B_tensor += B2_tensor;
    auto && S_voight = Material::stressToVoigt<dim>(S);
    r = B_tensor.transpose() * S_voight;
  }

  /// compute @f$ k_e = \int_e \mathbf{B}^t * \mathbf{D} * \mathbf{B}@f$
  auto r_e = std::make_unique<Array<Real>>(nb_element, bt_s_size, "r_e");
  fem.integrate(*bt_s, *r_e, bt_s_size, type, ghost_type, elem_filter);

  getModel().getDOFManager().assembleElementalArrayLocalArray(
      *r_e, internal_force, type, ghost_type, -1., elem_filter);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void Material::computePotentialEnergyByElements() {
  AKANTU_DEBUG_IN();

  for (auto type :
       getElementFilter().elementTypes(spatial_dimension, _not_ghost)) {
    computePotentialEnergy(type);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void Material::computePotentialEnergy(ElementType /*type*/) {
  AKANTU_TO_IMPLEMENT();
}

/* -------------------------------------------------------------------------- */
Real Material::getPotentialEnergy() {
  AKANTU_DEBUG_IN();
  Real epot = 0.;
  auto & fem = getFEEngine();
  computePotentialEnergyByElements();

  /// integrate the potential energy for each type of elements
  for (auto type :
       getElementFilter().elementTypes(spatial_dimension, _not_ghost)) {
    epot += fem.integrate(potential_energy(type, _not_ghost), type, _not_ghost,
                          getElementFilter(type, _not_ghost));
  }

  AKANTU_DEBUG_OUT();
  return epot;
}

/* -------------------------------------------------------------------------- */
Real Material::getPotentialEnergy(ElementType type, Int index) {
  return getPotentialEnergy({type, index, _not_ghost});
}

/* -------------------------------------------------------------------------- */
Real Material::getPotentialEnergy(const Element & element) {
  AKANTU_DEBUG_IN();
  auto & fem = getFEEngine();
  Vector<Real> epot_on_quad_points(fem.getNbIntegrationPoints(element.type));

  auto epot = fem.integrate(epot_on_quad_points,
                            {element.type,
                             getElementFilter(element.type)(element.element),
                             _not_ghost});

  AKANTU_DEBUG_OUT();
  return epot;
}

/* -------------------------------------------------------------------------- */
Real Material::getEnergy(const std::string & type) {
  AKANTU_DEBUG_IN();
  if (type == "potential") {
    return getPotentialEnergy();
  }
  AKANTU_DEBUG_OUT();
  return 0.;
}

/* -------------------------------------------------------------------------- */
Real Material::getEnergy(const std::string & energy_id,
                         const Element & element) {
  AKANTU_DEBUG_IN();
  if (energy_id == "potential") {
    return getPotentialEnergy(element);
  }
  AKANTU_DEBUG_OUT();
  return 0.;
}

/* -------------------------------------------------------------------------- */
void Material::initElementalFieldInterpolation(
    const ElementTypeMapArray<Real> & interpolation_points_coordinates) {
  AKANTU_DEBUG_IN();
  auto & fem = getFEEngine();
  fem.initElementalFieldInterpolationFromIntegrationPoints(
      interpolation_points_coordinates, this->interpolation_points_matrices,
      this->interpolation_inverse_coordinates, &(this->getElementFilter()));

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void Material::interpolateStress(ElementTypeMapArray<Real> & result,
                                 const GhostType ghost_type) {
  auto & fem = getFEEngine();
  fem.interpolateElementalFieldFromIntegrationPoints(
      this->stress, this->interpolation_points_matrices,
      this->interpolation_inverse_coordinates, result, ghost_type,
      &(this->getElementFilter()));
}

/* -------------------------------------------------------------------------- */
void Material::interpolateStressOnFacets(
    ElementTypeMapArray<Real> & result,
    ElementTypeMapArray<Real> & by_elem_result, const GhostType ghost_type) {
  interpolateStress(by_elem_result, ghost_type);

  auto stress_size = this->stress.getNbComponent();

  const auto & mesh = this->getModel().getMesh();
  const auto & mesh_facets = mesh.getMeshFacets();

  for (auto type :
       getElementFilter().elementTypes(spatial_dimension, ghost_type)) {
    auto nb_element_full = mesh.getNbElement(type, ghost_type);

    auto nb_interpolation_points_per_elem =
        by_elem_result(type, ghost_type).size() / nb_element_full;

    const auto & facet_to_element =
        mesh_facets.getSubelementToElement(type, ghost_type);

    auto nb_facet_per_elem = facet_to_element.getNbComponent();
    auto nb_quad_per_facet =
        nb_interpolation_points_per_elem / nb_facet_per_elem;
    Element element{type, 0, ghost_type};

    auto && by_elem_res =
        make_view(by_elem_result(type, ghost_type), stress_size,
                  nb_quad_per_facet, nb_facet_per_elem)
            .begin();

    for (auto global_el : getElementFilter(type, ghost_type)) {
      element.element = global_el;

      auto && result_per_elem = by_elem_res[global_el];

      for (Int f = 0; f < nb_facet_per_elem; ++f) {
        auto facet_elem = facet_to_element(global_el, f);
        Int is_second_element =
            Int(mesh_facets.getElementToSubelement(facet_elem)[0] != element);

        auto && result_local =
            result.get(facet_elem, stress_size, 2, nb_quad_per_facet);

        for (auto && data : zip(result_local, result_per_elem(f))) {
          std::get<0>(data)(is_second_element) = std::get<1>(data);
        }
      }
    }
  }
}

/* -------------------------------------------------------------------------- */
void Material::beforeSolveStep() { this->savePreviousState(); }

/* -------------------------------------------------------------------------- */
void Material::afterSolveStep(bool converged) {
  if (not converged) {
    this->restorePreviousState();
    return;
  }

  for (const auto & type : getElementFilter().elementTypes(
           _all_dimensions, _not_ghost, _ek_not_defined)) {
    this->updateEnergies(type);
  }
}
/* -------------------------------------------------------------------------- */
void Material::onDamageIteration() { this->savePreviousState(); }

/* -------------------------------------------------------------------------- */
void Material::onDamageUpdate() {
  for (const auto & type : getElementFilter().elementTypes(
           _all_dimensions, _not_ghost, _ek_not_defined)) {
    this->updateEnergiesAfterDamage(type);
  }
}

/* -------------------------------------------------------------------------- */
void Material::onDump() {
  if (this->isFiniteDeformation()) {
    this->computeAllCauchyStresses(_not_ghost);
  }
}

/* -------------------------------------------------------------------------- */
/// extrapolate internal values
void Material::extrapolateInternal(const ID & id, const Element & element,
                                   const Matrix<Real> & /*point*/,
                                   Matrix<Real> & extrapolated) {
  if (this->isInternal<Real>(id, element.kind())) {
    auto name = this->getID() + ":" + id;
    auto nb_quads =
        this->getInternal<Real>(name).getFEEngine().getNbIntegrationPoints(
            element.type, element.ghost_type);
    const auto & internal =
        this->getArray<Real>(id, element.type, element.ghost_type);
    auto nb_component = internal.getNbComponent();
    auto internal_it = make_view(internal, nb_component, nb_quads).begin();
    auto local_element = this->convertToLocalElement(element);

    /// instead of really extrapolating, here the value of the first GP
    /// is copied into the result vector. This works only for linear
    /// elements
    /// @todo extrapolate!!!!
    AKANTU_DEBUG_WARNING("This is a fix, values are not truly extrapolated");

    auto && values = internal_it[local_element.element];
    Idx index = 0;
    Vector<Real> tmp(nb_component);
    for (Int j = 0; j < values.cols(); ++j) {
      tmp = values(j);
      if (tmp.norm() > 0) {
        index = j;
        break;
      }
    }

    for (Int i = 0; i < extrapolated.size(); ++i) {
      extrapolated(i) = values(index);
    }
  } else {
    extrapolated.fill(0.);
  }
}

/* -------------------------------------------------------------------------- */
void Material::applyEigenGradU(const Matrix<Real> & prescribed_eigen_grad_u,
                               const GhostType ghost_type) {
  for (auto && type : getElementFilter().elementTypes(
           _all_dimensions, _not_ghost, _ek_not_defined)) {
    if (getElementFilter(type, ghost_type).empty()) {
      continue;
    }

    for (auto & eigengradu : make_view(this->eigengradu(type, ghost_type),
                                       spatial_dimension, spatial_dimension)) {
      eigengradu = prescribed_eigen_grad_u;
    }
  }
}

/* -------------------------------------------------------------------------- */
MaterialFactory & Material::getFactory() {
  return MaterialFactory::getInstance();
}

} // namespace akantu
