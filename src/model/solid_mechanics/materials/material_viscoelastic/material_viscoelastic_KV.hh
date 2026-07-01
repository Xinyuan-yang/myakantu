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

#ifndef AKANTU_MATERIAL_VISCOELASTIC_KELVIN_VOIGT_HH_
#define AKANTU_MATERIAL_VISCOELASTIC_KELVIN_VOIGT_HH_

namespace akantu {

/**
 * Material Viscoelastic based on Kelvin-Voigt model
 *
 *
 * @verbatim

              E_0
      ------|\/\/\|-------
      |                  |
   ---|                  |---
      |                  |
      -------------[|-----
                \Eta|

 @endverbatim
 *
 * keyword : viscoelastic_kelvin_voigt
 *
 * parameters in the material files :
 *   - Einf  : one spring element stiffness
 *   - Eta   : viscosity of the Kelvin-Voigt element
 */

template <Int dim>
class MaterialViscoelasticKelvinVoigt : public MaterialElastic<dim> {
  using Parent = MaterialElastic<dim>;
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  MaterialViscoelasticKelvinVoigt(SolidMechanicsModel & model, const ID & id = "");

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// initialize the material computed parameter
  void initMaterial() override;

  /// recompute the lame coefficient and effective tangent moduli
  void updateInternalParameters() override;

  /// update energy on a converged Newton
  void afterSolveStep(bool converged) override;
  
  /// constitutive law for all element of a type
  void computeStress(ElementType el_type,
                     GhostType ghost_type = _not_ghost) override;

  /// compute the tangent stiffness matrix for an element type
  void computeTangentModuli(ElementType el_type, Array<Real> & tangent_matrix,
                            GhostType ghost_type = _not_ghost) override;

  /// compute the elastic potential energy
  void computePotentialEnergy(ElementType el_type) override;

protected:
  template <typename D>
  void computePotentialEnergyOnQuad(const Eigen::MatrixBase<D> & grad_u,
                                    Real & epot);

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
    return MaterialElastic<dim>::getArguments(el_type, ghost_type);
  }

  MatrixType getTangentType() const override { return _symmetric; }

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

  /// Viscosity of the Kelvin-Voigt dashpot
  Real Eta{0.};

  Real Einf{0.};

  /// time step from previous solveStep
  Real previous_dt{0.};

  /// Stiffness matrix template
  Matrix<Real, voigt_h::size, voigt_h::size> C;
  /// Compliance matrix template
  Matrix<Real, voigt_h::size, voigt_h::size> D;

  /// Dissipated energy
  InternalField<Real> & dissipated_energy;

  /// Mechanical work
  InternalField<Real> & mechanical_work;

};

} // namespace akantu

#endif /* AKANTU_MATERIAL_VISCOELASTIC_KELVIN_VOIGT_HH_ */
