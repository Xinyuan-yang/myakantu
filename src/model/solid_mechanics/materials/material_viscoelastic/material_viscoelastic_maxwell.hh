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
#include "aka_common.hh"
#include "aka_voigthelper.hh"
#include "material_elastic.hh"
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_MATERIAL_VISCOELASTIC_MAXWELL_HH_
#define AKANTU_MATERIAL_VISCOELASTIC_MAXWELL_HH_

namespace akantu {

/**
 * Material Viscoelastic based on Maxwell chain
 *
 *
 * @verbatim

              E_0
      ------|\/\/\|-------
      |                  |
   ---|                  |---
      |                  |
      ----|\/\/\|--[|-----
      |   E_v1  \Eta_1|
   ---|                  |---
      |                  |
      ----|\/\/\|--[|-----
      |   E_v2 \Eta_2 |
   ---|                  |---
      |                  |
      ----|\/\/\|--[|----
          E_vN \Eta_N

 @endverbatim
 *
 * keyword : viscoelastic_maxwell
 *
 * parameters in the material files :
 *   - N   : number of Maxwell elements
 *   - Einf  : one spring element stiffness
 *   - Ev1 : stiffness of the 1st viscous element
 *   - Eta1: viscosity of the 1st Maxwell element
 *   ...
 *   - Ev<N> : stiffness of the Nst viscous element
 *   - Eta<N>: viscosity of the Nst Maxwell element
 */

template <Int dim>
class MaterialViscoelasticMaxwell : public MaterialElastic<dim> {
  using Parent = MaterialElastic<dim>;
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  MaterialViscoelasticMaxwell(SolidMechanicsModel & model, const ID & id = "");

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// initialize the material computed parameter
  void initMaterial() override;

  /// recompute the lame coefficient and effective tangent moduli
  void updateInternalParameters() override;

  /// update internal variable on a converged Newton
  void afterSolveStep(bool converged) override;

  /// update internal variable based on previous and current strain values
  void updateIntVariables();

  /// update the internal variable sigma_v on quadrature point
  template <typename D1, typename D2>
  void updateIntVarOnQuad(const Eigen::MatrixBase<D1> & grad_u,
                          const Eigen::MatrixBase<D2> & previous_grad_u,
                          Tensor3Proxy<Real> & sigma_v,
                          Tensor3Proxy<Real> & epsilon_v);

  /// constitutive law for all element of a type
  void computeStress(ElementType el_type,
                     GhostType ghost_type = _not_ghost) override;

  /// compute the tangent stiffness matrix for an element type
  void computeTangentModuli(ElementType el_type, Array<Real> & tangent_matrix,
                            GhostType ghost_type = _not_ghost) override;

  /// change flag of updateIntVar to true
  void forceUpdateVariable();

  /// change flag of updateIntVar to false
  void forceNotUpdateVariable();

  /// compute the elastic potential energy
  void computePotentialEnergy(ElementType el_type) override;

protected:
  template <typename D>
  void computePotentialEnergyOnQuad(const Eigen::MatrixBase<D> & grad_u,
                                    Real & epot, Tensor3Proxy<Real> & sigma_v,
                                    Tensor3Proxy<Real> & epsilon_v);

  /// update the dissipated energy, is called after the stress have been
  /// computed
  void updateDissipatedEnergy(ElementType el_type);

  template <typename D1, typename D2, typename D3, typename D4>
  void
  updateDissipatedEnergyOnQuad(const Eigen::MatrixBase<D1> & grad_u,
                               const Eigen::MatrixBase<D2> & previous_grad_u,
                               const Eigen::MatrixBase<D3> & sigma,
                               const Eigen::MatrixBase<D4> & previous_sigma,
                               Real & dis_energy, Real & mech_work,
                               const Real & pot_energy);

  /// compute stresses on a quadrature point
  template <class Args> void computeStressOnQuad(Args && args);

  /// compute tangent moduli on a quadrature point
  template <typename D1>
  void computeTangentModuliOnQuad(Eigen::MatrixBase<D1> & tangent);

  bool hasStiffnessMatrixChanged() const override {
    Real dt = this->getModel().getTimeStep();

    return ((this->previous_dt == dt)
                ? (!(this->previous_dt == dt)) * (this->was_stiffness_assembled)
                : (!(this->previous_dt == dt)));
    //  return (!(this->previous_dt == dt));
  }

  decltype(auto) getArguments(ElementType el_type,
                              GhostType ghost_type = _not_ghost) {
    return zip_append(MaterialElastic<dim>::getArguments(el_type, ghost_type),
                      "sigma_v"_n = make_view((*sigma_v)(el_type, ghost_type),
                                              dim, dim, Ev.size()));
  }

  MatrixType getTangentType() override { return _symmetric; }

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  /// give the dissipated energy
  [[nodiscard]] Real getDissipatedEnergy() const;
  [[nodiscard]] Real getDissipatedEnergy(const Element & element) const;

  /// get the potential energy
  [[nodiscard]] Real getPotentialEnergy() const;
  [[nodiscard]] Real getPotentialEnergy(const Element & element) const;

  /// get the potential energy
  [[nodiscard]] Real getMechanicalWork() const;
  [[nodiscard]] Real getMechanicalWork(const Element & element) const;

  /// get the energy using an energy type string for the time step
  [[nodiscard]] Real getEnergy(const std::string & type) override;
  [[nodiscard]] Real getEnergy(const std::string & energy_id,
                               const Element & element) override;
  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  using voigt_h = VoigtHelper<dim>;

  /// Vectors of viscosity, viscous elastic modulus, one spring element elastic
  /// modulus
  Vector<Real> Eta;
  Vector<Real> Ev;

  Real Einf{0.};

  /// time step from previous solveStep
  Real previous_dt{0.};

  /// Stiffness matrix template
  Matrix<Real, voigt_h::size, voigt_h::size> C;
  /// Compliance matrix template
  Matrix<Real, voigt_h::size, voigt_h::size> D;

  /// Internal variable: viscous_stress
  std::shared_ptr<InternalField<Real>> sigma_v;

  /// Internal variable: spring strain in Maxwell element
  std::shared_ptr<InternalField<Real>> epsilon_v;

  /// Dissipated energy
  InternalField<Real> & dissipated_energy;

  /// Mechanical work
  InternalField<Real> & mechanical_work;

  /// Update internal variable after solve step or not
  bool update_variable_flag{true};
};

} // namespace akantu

#endif /* AKANTU_MATERIAL_VISCOELASTIC_MAXWELL_HH_ */
