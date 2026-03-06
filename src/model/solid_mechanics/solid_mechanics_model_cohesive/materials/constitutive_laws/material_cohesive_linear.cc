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
#include "material_cohesive_linear.hh"
#include "aka_array.hh"
#include "aka_common.hh"
#include "dof_synchronizer.hh"
#include "solid_mechanics_model_cohesive.hh"
#include "sparse_matrix.hh"
/* -------------------------------------------------------------------------- */
#include <algorithm>
#include <cmath>
#include <vector>
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
template <Int dim>
MaterialCohesiveLinear<dim>::MaterialCohesiveLinear(SolidMechanicsModel & model,
                                                    const ID & id)
    : MaterialCohesive(model, id),
      sigma_c_eff(registerInternal<Real, CohesiveRandomInternalField>(
          "sigma_c_eff", 1)),
      delta_c_eff(
          registerInternal<Real, CohesiveInternalField>("delta_c_eff", 1)),
      insertion_stress(registerInternal<Real, CohesiveInternalField>(
          "insertion_stress", dim)) {
  AKANTU_DEBUG_IN();

  this->registerParam("beta", beta, Real(0.), _pat_parsable | _pat_readable,
                      "Beta parameter");

  this->registerParam("G_c", G_c, Real(0.), _pat_parsable | _pat_readable,
                      "Mode I fracture energy");

  this->registerParam("penalty", penalty, Real(0.),
                      _pat_parsable | _pat_readable, "Penalty coefficient");

  this->registerParam("volume_s", volume_s, Real(0.),
                      _pat_parsable | _pat_readable,
                      "Reference volume for sigma_c scaling");

  this->registerParam("m_s", m_s, Real(1.), _pat_parsable | _pat_readable,
                      "Weibull exponent for sigma_c scaling");

  this->registerParam("kappa", kappa, Real(1.), _pat_parsable | _pat_readable,
                      "Kappa parameter");

  this->registerParam(
      "contact_after_breaking", contact_after_breaking, false,
      _pat_parsable | _pat_readable,
      "Activation of contact when the elements are fully damaged");

  this->registerParam("max_quad_stress_insertion", max_quad_stress_insertion,
                      false, _pat_parsable | _pat_readable,
                      "Insertion of cohesive element when stress is high "
                      "enough just on one quadrature point");

  this->registerParam("recompute", recompute, false,
                      _pat_parsable | _pat_modifiable, "recompute solution");

  this->use_previous_delta_max = true;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <Int dim> void MaterialCohesiveLinear<dim>::initMaterial() {
  AKANTU_DEBUG_IN();

  MaterialCohesive::initMaterial();

  if (not Math::are_float_equal(delta_c, 0.)) {
    delta_c_eff.setDefaultValue(delta_c);
  } else {
    Real sigma_c = this->sigma_c;
    delta_c_eff.setDefaultValue(2 * G_c / sigma_c);
  }

  if (model->getIsExtrinsic()) {
    scaleInsertionTraction();
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <Int dim>
void MaterialCohesiveLinear<dim>::updateInternalParameters() {
  /// compute scalars
  beta2_kappa2 = beta * beta / kappa / kappa;
  beta2_kappa = beta * beta / kappa;

  if (Math::are_float_equal(beta, 0)) {
    beta2_inv = 0;
  } else {
    beta2_inv = 1. / beta / beta;
  }
}

/* -------------------------------------------------------------------------- */
template <Int dim> void MaterialCohesiveLinear<dim>::scaleInsertionTraction() {
  AKANTU_DEBUG_IN();

  // do nothing if volume_s hasn't been specified by the user
  if (Math::are_float_equal(volume_s, 0.)) {
    return;
  }

  const auto & mesh_facets = model->getMeshFacets();
  const auto & fe_engine = model->getFEEngine();
  const auto & fe_engine_facet = model->getFEEngine("FacetsFEEngine");

  for (const auto & type_facet : mesh_facets.elementTypes(dim - 1)) {
    const auto & facet_to_element =
        mesh_facets.getElementToSubelement(type_facet);

    auto nb_quad_per_facet = fe_engine_facet.getNbIntegrationPoints(type_facet);

    for (auto data :
         enumerate(make_view(sigma_c(type_facet), nb_quad_per_facet))) {
      auto f = std::get<0>(data);
      auto && sigma_c = std::get<1>(data);

      // compute bounding volume
      Real volume = 0;

      for (auto && elem : facet_to_element(f)) {
        if (elem == ElementNull) {
          continue;
        }

        // unit vector for integration in order to obtain the volume
        auto nb_quadrature_points = fe_engine.getNbIntegrationPoints(elem.type);
        Vector<Real> unit_vector(nb_quadrature_points);
        unit_vector.fill(1);

        volume += fe_engine.integrate(unit_vector, elem);
      }

      // scale sigma_c
      Vector<Real> base_sigma_c_v(sigma_c.rows());
      sigma_c = (sigma_c.colwise() - base_sigma_c_v) *
                    std::pow(volume_s / volume, 1. / m_s) +
                base_sigma_c_v;
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <Int dim> // NOLINTNEXTLINE(readability-function-cognitive-complexity)
void MaterialCohesiveLinear<dim>::checkInsertion(bool check_only) {
  AKANTU_DEBUG_IN();

  const auto & mesh_facets = model->getMeshFacets();
  auto & inserter = model->getElementInserter();

  for (const auto & type_facet : mesh_facets.elementTypes(dim - 1)) {
    auto type_cohesive = FEEngine::getCohesiveElementType(type_facet);

    const auto & facets_check = inserter.getCheckFacets(type_facet);
    auto & f_insertion = inserter.getInsertionFacets(type_facet);
    auto & sig_c_eff = sigma_c_eff(type_cohesive);
    auto & del_c = delta_c_eff(type_cohesive);
    auto & ins_stress = insertion_stress(type_cohesive);
    auto & trac_old = tractions.previous(type_cohesive);
    const auto & f_stress = model->getStressOnFacets(type_facet);

    const auto & facet_filter_array = getFacetFilter(type_facet);
    const auto & sigma_limit_array = sigma_c(type_facet);

    auto nb_quad_facet =
        model->getFEEngine("FacetsFEEngine").getNbIntegrationPoints(type_facet);

#ifndef AKANTU_NDEBUG
    auto nb_quad_cohesive = model->getFEEngine("CohesiveFEEngine")
                                .getNbIntegrationPoints(type_cohesive);

    AKANTU_DEBUG_ASSERT(nb_quad_cohesive == nb_quad_facet,
                        "The cohesive element and the corresponding facet do "
                        "not have the same numbers of integration points");
#endif

    Matrix<Real> normal_traction(dim, nb_quad_facet);
    Vector<Real> stress_check(nb_quad_facet);

    const auto & tangents = model->getTangents(type_facet);
    const auto & normals = model->getFEEngine("FacetsFEEngine")
                               .getNormalsOnIntegrationPoints(type_facet);
    auto normal_begin = make_view<dim>(normals).begin();
    auto tangent_begin = make_view<dim, (dim == 3 ? 2 : 1)>(tangents).begin();
    auto facet_stress_begin = make_view(f_stress, dim, dim, 2).begin();

    std::vector<Real> new_sigmas;
    std::vector<Vector<Real>> new_normal_traction;
    std::vector<Real> new_delta_c;

    // loop over each facet belonging to this material
    for (auto && [facet, sigma_limit] :
         zip(facet_filter_array, sigma_limit_array)) {
      // skip facets where check shouldn't be realized
      if (!facets_check(facet)) {
        continue;
      }

      // compute the effective norm on each quadrature point of the facet
      for (Int q = 0; q < nb_quad_facet; ++q) {
        auto current_quad = (facet * nb_quad_facet) + q;
        auto && normal = normal_begin[current_quad];
        auto && tangent = tangent_begin[current_quad];
        auto && facet_stress = facet_stress_begin[current_quad];

        // compute average stress on the current quadrature point
        auto && stress_1 = facet_stress(0);
        auto && stress_2 = facet_stress(1);

        auto && stress_avg = (stress_1 + stress_2) / 2.;

        // compute normal and effective stress
        stress_check(q) = computeEffectiveNorm(stress_avg, normal, tangent,
                                               normal_traction(q));
      }

      // verify if the effective stress overcomes the threshold
      auto final_stress = stress_check.mean();
      if (max_quad_stress_insertion) {
        final_stress =
            *std::max_element(stress_check.begin(), stress_check.end());
      }

      if (final_stress > sigma_limit) {
        f_insertion(facet) = true;

        if (check_only) {
          continue;
        }

        // store the new cohesive material parameters for each quadrature
        // point
        for (Int q = 0; q < nb_quad_facet; ++q) {
          auto new_sigma = stress_check(q);
          auto && normal_traction_vec = normal_traction(q);

          if (dim != 3) {
            normal_traction_vec *= -1.;
          }

          new_sigmas.push_back(new_sigma);
          new_normal_traction.emplace_back(normal_traction_vec);

          Real new_delta{};

          // set delta_c in function of G_c or a given delta_c value
          if (Math::are_float_equal(delta_c, 0.)) {
            new_delta = 2 * G_c / new_sigma;
          } else {
            new_delta = sigma_limit / new_sigma * delta_c;
          }

          new_delta_c.push_back(new_delta);
        }
      }
    }

    // update material data for the new elements
    auto old_nb_quad_points = sig_c_eff.size();
    Int new_nb_quad_points = Int(new_sigmas.size());
    sig_c_eff.resize(old_nb_quad_points + new_nb_quad_points);
    ins_stress.resize(old_nb_quad_points + new_nb_quad_points);
    trac_old.resize(old_nb_quad_points + new_nb_quad_points);
    del_c.resize(old_nb_quad_points + new_nb_quad_points);

    for (Int q = 0; q < new_nb_quad_points; ++q) {
      sig_c_eff(old_nb_quad_points + q) = new_sigmas[q];
      del_c(old_nb_quad_points + q) = new_delta_c[q];
      for (Int d = 0; d < dim; ++d) {
        ins_stress(old_nb_quad_points + q, d) = new_normal_traction[q](d);
        trac_old(old_nb_quad_points + q, d) = new_normal_traction[q](d);
      }
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <Int dim>
void MaterialCohesiveLinear<dim>::computeTraction(ElementType el_type,
                                                  GhostType ghost_type) {
  for (auto && args : getArguments(el_type, ghost_type)) {
    this->computeTractionOnQuad(args);
  }
}

/* -------------------------------------------------------------------------- */
template <Int dim>
void MaterialCohesiveLinear<dim>::computeTangentTraction(
    ElementType el_type, Array<Real> & tangent_matrix, GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  for (auto && [args, tangent] : zip(getArguments(el_type, ghost_type),
                                     make_view<dim, dim>(tangent_matrix))) {
    computeTangentTractionOnQuad(tangent, args);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template class MaterialCohesiveLinear<1>;
template class MaterialCohesiveLinear<2>;
template class MaterialCohesiveLinear<3>;
const bool material_is_alocated_cohesive_linear [[maybe_unused]] =
    instantiateMaterial<MaterialCohesiveLinear>("cohesive_linear");

} // namespace akantu
