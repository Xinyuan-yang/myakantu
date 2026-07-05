/**
 * Copyright (©) 2018-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "material_viscoelastic_KV.hh"
#include "solid_mechanics_model.hh"

namespace akantu {

/* -------------------------------------------------------------------------- */
template <Int dim>
MaterialViscoelasticKelvinVoigt<dim>::MaterialViscoelasticKelvinVoigt(
    SolidMechanicsModel & model, const ID & id)
    : MaterialElastic<dim>(model, id),
      dissipated_energy(this->registerInternal("dissipated_energy", 1)),
      mechanical_work(this->registerInternal("mechanical_work", 1)) {
  this->registerParam("Einf", Einf, Real(1.), _pat_parsable | _pat_modifiable,
                      "Stiffness of the elastic element");
  this->registerParam("previous_dt", previous_dt, Real(0.), _pat_readable,
                      "Time step of previous solveStep");
  this->registerParam("Eta", Eta, Real(0.), _pat_parsable | _pat_modifiable,
                      "Viscosity of element");
}

/* -------------------------------------------------------------------------- */
template <Int dim> void MaterialViscoelasticKelvinVoigt<dim>::initMaterial() {
  AKANTU_DEBUG_IN();

  this->E = Einf;

  MaterialElastic<dim>::initMaterial();

  AKANTU_DEBUG_ASSERT(
      !this->finite_deformation,
      "Current material works only in infinitesimal deformations.");

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <Int dim>
void MaterialViscoelasticKelvinVoigt<dim>::updateInternalParameters() {
  MaterialElastic<dim>::updateInternalParameters();

  [[maybe_unused]] auto pre_mult = 1 / (1 + this->nu) / (1 - 2 * this->nu);

  [[maybe_unused]] auto Miiii = pre_mult * (1 - this->nu);
  [[maybe_unused]] auto Miijj = pre_mult * this->nu;
  [[maybe_unused]] auto Mijij = pre_mult * 0.5 * (1 - 2 * this->nu);

  if constexpr (dim == 2) {
    if (this->plane_stress) {
      auto plane_stress_pre_mult = 1 / (1 - this->nu * this->nu);
      Miiii = plane_stress_pre_mult;
      Miijj = plane_stress_pre_mult * this->nu;
      Mijij = plane_stress_pre_mult * 0.5 * (1 - this->nu);
    }
  }

  [[maybe_unused]] auto Diiii = 1;
  [[maybe_unused]] auto Diijj = -this->nu;
  [[maybe_unused]] auto Dijij = (2 + 2 * this->nu);

  C.zero();
  D.zero();

  if constexpr (dim == 1) {
    C(0, 0) = 1;
    D(0, 0) = 1;
  } else {
    C(0, 0) = Miiii;
    D(0, 0) = Diiii;
  }
  if constexpr (dim >= 2) {
    auto n = voigt_h::size;
    C(1, 1) = Miiii;
    C(0, 1) = Miijj;
    C(1, 0) = Miijj;
    C(n - 1, n - 1) = Mijij;

    D(1, 1) = Diiii;
    D(0, 1) = Diijj;
    D(1, 0) = Diijj;
    D(n - 1, n - 1) = Dijij;
  }

  if constexpr (dim == 3) {
    C(2, 2) = Miiii;
    C(0, 2) = Miijj;
    C(1, 2) = Miijj;
    C(2, 0) = Miijj;
    C(2, 1) = Miijj;
    C(3, 3) = Mijij;
    C(4, 4) = Mijij;

    D(2, 2) = Diiii;
    D(0, 2) = Diijj;
    D(1, 2) = Diijj;
    D(2, 0) = Diijj;
    D(2, 1) = Diijj;
    D(3, 3) = Dijij;
    D(4, 4) = Dijij;
  }
}

/* -------------------------------------------------------------------------- */
template <Int dim>
void MaterialViscoelasticKelvinVoigt<dim>::computeStress(ElementType el_type,
                                                     GhostType ghost_type) {
  // NOLINTNEXTLINE(bugprone-parent-virtual-call)
  Parent::computeStress(el_type, ghost_type);

  for (auto && args : getArguments(el_type, ghost_type)) {
    computeStressOnQuad(args);
  }
}

/* -------------------------------------------------------------------------- */
template <Int dim>
template <class Args>
void MaterialViscoelasticKelvinVoigt<dim>::computeStressOnQuad(Args && args) {
  Real dt = this->getModel().getTimeStep();
  // Wikipedia convention:
  // 2*eps_ij (i!=j) = voigt_eps_I
  // http://en.wikipedia.org/wiki/Voigt_notation
  auto voigt_current_strain = Material::strainToVoigt<dim>(
      Material::gradUToEpsilon<dim>(args["grad_u"_n]));
  auto voigt_previous_strain = Material::strainToVoigt<dim>(
      Material::gradUToEpsilon<dim>(args["previous_grad_u"_n]));

  Vector<Real, voigt_h::size> voigt_stress =
      this->Einf * this->C * voigt_current_strain;

  Vector<Real, voigt_h::size> voigt_strain_increment =
      this->C * (voigt_current_strain - voigt_previous_strain);

  voigt_stress += this->Eta / dt * voigt_strain_increment;

  auto && sigma = args["sigma"_n];
  for (Int I = 0; I < voigt_h::size; ++I) {
    auto && [i, j] = voigt_h::vec[I];
    sigma(i, j) = sigma(j, i) =
        voigt_stress(I) + Math::kronecker(i, j) * args["sigma_th"_n];
  }
}

/* -------------------------------------------------------------------------- */
template <Int dim>
void MaterialViscoelasticKelvinVoigt<dim>::computePotentialEnergy(
    ElementType el_type) {
  for (auto && [args, epot] :
       zip(getArguments(el_type), this->potential_energy(el_type))) {
    this->computePotentialEnergyOnQuad(args["grad_u"_n], epot);
  }
}

/* -------------------------------------------------------------------------- */
template <Int dim>
template <typename D1>
void MaterialViscoelasticKelvinVoigt<dim>::computePotentialEnergyOnQuad(
    const Eigen::MatrixBase<D1> & grad_u, Real & epot) {
  auto voigt_strain =
      Material::strainToVoigt<dim>(Material::gradUToEpsilon<dim>(grad_u));
  Vector<Real, voigt_h::size> voigt_stress =
      this->Einf * this->C * voigt_strain;

  epot = 0.5 * voigt_stress.dot(voigt_strain);
}

/* -------------------------------------------------------------------------- */
template <Int dim>
void MaterialViscoelasticKelvinVoigt<dim>::afterSolveStep(bool converged) {
  Material::afterSolveStep(converged);

  if (not converged) {
    return;
  }

  for (const auto & el_type : this->getElementFilter().elementTypes(
           _all_dimensions, _not_ghost, _ek_not_defined)) {
    this->updateDissipatedEnergy(el_type);
  }
}
/* -------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------- */
template <Int dim>
void MaterialViscoelasticKelvinVoigt<dim>::computeTangentModuli(
    ElementType /*el_type*/, Array<Real> & tangent_matrix,
    GhostType /*ghost_type*/) {
  AKANTU_DEBUG_IN();

  Real dt = this->getModel().getTimeStep();
  Real E_ef = this->Einf + this->Eta / dt;

  this->previous_dt = dt;

  const auto tangent_size = Material::getTangentStiffnessVoigtSize(dim);
  for (auto && tangent :
       make_view<tangent_size, tangent_size>(tangent_matrix)) {
    this->computeTangentModuliOnQuad(tangent);
  }

  tangent_matrix *= E_ef;

  this->was_stiffness_assembled = true;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <Int dim>
template <typename D1>
void MaterialViscoelasticKelvinVoigt<dim>::computeTangentModuliOnQuad(
    Eigen::MatrixBase<D1> & tangent) {

  tangent = C;
}

template <Int dim>
void MaterialViscoelasticKelvinVoigt<dim>::updateDissipatedEnergy(
    ElementType el_type) {
  this->computePotentialEnergy(el_type);

  for (auto && [args, epot, edis, work] :
       zip(getArguments(el_type), this->potential_energy(el_type),
           this->dissipated_energy(el_type), this->mechanical_work(el_type))) {
    updateDissipatedEnergyOnQuad(args["grad_u"_n], args["previous_grad_u"_n],
                                 args["sigma"_n], args["previous_sigma"_n],
                                 edis, work, epot);
  }
}

/* -------------------------------------------------------------------------- */
template <Int dim>
template <typename D1, typename D2, typename D3, typename D4>
void MaterialViscoelasticKelvinVoigt<dim>::updateDissipatedEnergyOnQuad(
    const Eigen::MatrixBase<D1> & grad_u,
    const Eigen::MatrixBase<D2> & previous_grad_u,
    const Eigen::MatrixBase<D3> & sigma,
    const Eigen::MatrixBase<D4> & previous_sigma, Real & dis_energy,
    Real & mech_work, const Real & pot_energy) {
  Real dt = this->getModel().getTimeStep();

  auto && strain_rate = (grad_u - previous_grad_u) / dt;
  auto && av_stress = (sigma + previous_sigma) / 2.;

  mech_work += av_stress.doubleDot(strain_rate) * dt;

  dis_energy = mech_work - pot_energy;
}

/* -------------------------------------------------------------------------- */
template <Int dim>
Real MaterialViscoelasticKelvinVoigt<dim>::getDissipatedEnergy() const {
  AKANTU_DEBUG_IN();

  auto & fem = this->getFEEngine();
  Real de = 0.;

  /// integrate the dissipated energy for each type of elements
  for (auto && type : this->getElementFilter().elementTypes(dim, _not_ghost)) {
    de += fem.integrate(this->dissipated_energy(type, _not_ghost), type,
                        _not_ghost, this->getElementFilter(type, _not_ghost));
  }

  AKANTU_DEBUG_OUT();
  return de;
}

/* -------------------------------------------------------------------------- */
template <Int dim>
Real MaterialViscoelasticKelvinVoigt<dim>::getDissipatedEnergy(
    const Element & element) const {
  auto & fem = this->getFEEngine();
  auto nb_quadrature_points = fem.getNbIntegrationPoints(element.type);
  auto it =
      make_view(this->dissipated_energy(element.type), nb_quadrature_points)
          .begin();
  auto mat_element = element;
  mat_element.element = this->getElementFilter()(element);

  return fem.integrate(it[element.element], mat_element);
}

/* -------------------------------------------------------------------------- */
template <Int dim>
Real MaterialViscoelasticKelvinVoigt<dim>::getMechanicalWork() const {
  auto & fem = this->getFEEngine();
  Real mw = 0.;

  /// integrate the dissipated energy for each type of elements
  for (auto && type : this->getElementFilter().elementTypes(dim)) {
    mw += fem.integrate(this->mechanical_work(type), type, _not_ghost,
                        this->getElementFilter(type, _not_ghost));
  }

  return mw;
}

/* -------------------------------------------------------------------------- */
template <Int dim>
Real MaterialViscoelasticKelvinVoigt<dim>::getMechanicalWork(
    const Element & element) const {
  auto & fem = this->getFEEngine();
  auto nb_quadrature_points = fem.getNbIntegrationPoints(element.type);
  auto it = make_view(this->mechanical_work(element.type), nb_quadrature_points)
                .begin();

  auto mat_element = element;
  mat_element.element = this->getElementFilter()(element);

  return fem.integrate(it[element.element], mat_element);
}

/* -------------------------------------------------------------------------- */
template <Int dim>
Real MaterialViscoelasticKelvinVoigt<dim>::getPotentialEnergy() const {
  AKANTU_DEBUG_IN();

  auto & fem = this->getFEEngine();
  Real epot = 0.;

  /// integrate the dissipated energy for each type of elements
  for (auto && type : this->getElementFilter().elementTypes(dim, _not_ghost)) {
    epot += fem.integrate(this->potential_energy(type, _not_ghost), type,
                          _not_ghost, this->getElementFilter(type, _not_ghost));
  }

  AKANTU_DEBUG_OUT();
  return epot;
}

/* -------------------------------------------------------------------------- */
template <Int dim>
Real MaterialViscoelasticKelvinVoigt<dim>::getPotentialEnergy(
    const Element & element) const {
  AKANTU_DEBUG_IN();

  auto & fem = this->getFEEngine();
  auto nb_quadrature_points = fem.getNbIntegrationPoints(element.type);
  auto it =
      make_view(this->potential_energy(element.type), nb_quadrature_points)
          .begin();
  auto mat_element = element;
  mat_element.element = this->getElementFilter()(element);

  AKANTU_DEBUG_OUT();
  return fem.integrate(it[element.element], mat_element);
}

/* -------------------------------------------------------------------------- */
template <Int dim>
Real MaterialViscoelasticKelvinVoigt<dim>::getEnergy(const std::string & type) {
  if (type == "dissipated") {
    return getDissipatedEnergy();
  }
  if (type == "potential") {
    return getPotentialEnergy();
  }
  if (type == "work") {
    return getMechanicalWork();
  }
  return MaterialElastic<dim>::getEnergy(type);
}

/* -------------------------------------------------------------------------- */
template <Int dim>
Real MaterialViscoelasticKelvinVoigt<dim>::getEnergy(const std::string & energy_id,
                                                 const Element & element) {
  if (energy_id == "dissipated") {
    return getDissipatedEnergy(element);
  }
  if (energy_id == "potential") {
    return getPotentialEnergy(element);
  }
  if (energy_id == "work") {
    return getMechanicalWork(element);
  }
  return MaterialElastic<dim>::getEnergy(energy_id, element);
}

/* -------------------------------------------------------------------------- */
const bool material_is_allocated_viscoelastic_kelvin_voigt [[maybe_unused]] =
    instantiateMaterial<MaterialViscoelasticKelvinVoigt>("viscoelastic_kelvin_voigt");

} // namespace akantu
