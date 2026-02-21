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
#include "plane_stress_toolbox.hh"

/* -------------------------------------------------------------------------- */

#ifndef AKANTU_MATERIAL_NEOHOOKEAN_HH_
#define AKANTU_MATERIAL_NEOHOOKEAN_HH_

namespace akantu {

/**
 * Material elastic isotropic
 *
 * parameters in the material files :
 *   - rho : density (default: 0)
 *   - E   : Young's modulus (default: 0)
 *   - nu  : Poisson's ratio (default: 1/2)
 *   - Plane_Stress : if 0: plane strain, else: plane stress (default: 0)
 */
template <Int dim> class MaterialNeohookean : public PlaneStressToolbox<dim> {
  using Parent = PlaneStressToolbox<dim>;
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  MaterialNeohookean(SolidMechanicsModel & model, const ID & id = "");

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// initialize the material computed parameter
  void initMaterial() override;

  /// constitutive law for all element of a type
  void computeStress(ElementType el_type,
                     GhostType ghost_type = _not_ghost) override;

  /// Computation of the cauchy stress for plane strain materials
  void
  computeCauchyStressPlaneStress(ElementType el_type,
                                 GhostType ghost_type = _not_ghost) override;

  /// Non linear computation of the third direction strain in 2D plane stress
  /// case
  void computeThirdAxisDeformation(ElementType el_type,
                                   GhostType ghost_type = _not_ghost) override;

  /// compute the elastic potential energy
  void computePotentialEnergy(ElementType el_type) override;

  /// compute the tangent stiffness matrix for an element type
  void computeTangentModuli(ElementType el_type, Array<Real> & tangent_matrix,
                            GhostType ghost_type = _not_ghost) override;

  /// compute the p-wave speed in the material
  [[nodiscard]] Real getPushWaveSpeed(const Element & element) const override;

  /// compute the s-wave speed in the material
  [[nodiscard]] Real getShearWaveSpeed(const Element & element) const override;

  MatrixType getTangentType() const override { return _symmetric; }

protected:
  /// constitutive law for a given quadrature point
  inline void computePiolaKirchhoffOnQuad(const Matrix<Real> & E,
                                          Matrix<Real> & S);

  /// constitutive law for a given quadrature point (first piola)
  inline void computeFirstPiolaKirchhoffOnQuad(const Matrix<Real> & grad_u,
                                               const Matrix<Real> & S,
                                               Matrix<Real> & P);

  /// constitutive law for a given quadrature point
  inline void computeDeltaStressOnQuad(const Matrix<Real> & grad_u,
                                       const Matrix<Real> & grad_delta_u,
                                       Matrix<Real> & delta_S);

  /// constitutive law for a given quadrature point
  template <class Args> inline void computeStressOnQuad(Args && args);

  /// constitutive law for a given quadrature point
  template <class Args>
  inline void computeThirdAxisDeformationOnQuad(Args && args);

  /// constitutive law for a given quadrature point
  // inline void updateStressOnQuad(const Matrix<Real> & sigma,
  //                              Matrix<Real> & cauchy_sigma);

  /// compute the potential energy for a quadrature point
  template <class Args>
  inline void computePotentialEnergyOnQuad(Args && args, Real & epot);

  /// compute the tangent stiffness matrix for an element
  template <class Args> void computeTangentModuliOnQuad(Args && args);

  /// recompute the lame coefficient if E or nu changes
  void updateInternalParameters() override;

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  /// the young modulus
  Real E{0.};

  /// Poisson coefficient
  Real nu{0.};

  /// First Lamé coefficient
  Real lambda{0.};

  /// Second Lamé coefficient (shear modulus)
  Real mu{0.};

  /// Bulk modulus
  Real kpa{0.};
};

} // namespace akantu

/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */
#include "material_neohookean_inline_impl.hh"

#endif /* AKANTU_MATERIAL_NEOHOOKEAN_HH_ */
