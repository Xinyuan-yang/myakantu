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
#ifndef __AST_NTN_FRICLAW_RATE_AND_STATE_DIETERICH_HH__
#define __AST_NTN_FRICLAW_RATE_AND_STATE_DIETERICH_HH__

/* -------------------------------------------------------------------------- */
// simtools
#include "ntn_friclaw_coulomb.hh"

namespace akantu{
  
/* -------------------------------------------------------------------------- */
template <class Regularisation = NTNFricRegNoRegularisation>
class NTNFricLawRateAndStateDieterich : public NTNFricLawCoulomb<Regularisation> {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  
  NTNFricLawRateAndStateDieterich(NTNBaseContact & contact,
				  const ID & id = "rate_and_state_dieterich");
  virtual ~NTNFricLawRateAndStateDieterich() {};
  
  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// register synchronizedarrays for sync
  void registerSynchronizedArray(SynchronizedArrayBase & array) override;

  /// dump restart file
  void dumpRestart(const std::string & file_name) const override;

  /// read restart file
  void readRestart(const std::string & file_name) override;

  /// function to print the contain of the class
  void printself(std::ostream & stream, int indent = 0) const override;
  /// add noise to the state variable. A map is passed to each process containing the full noise profil to avoid any miscomunication errors
  void addNoiseToStateVariable(std::map<Real,Real> &noise);
  
protected:
  /// compute frictional strength according to friction law
  void computeFrictionalStrength() override;
  /// computes the friction coefficient as a function of slip
  void computeFrictionCoefficient();

  /* ------------------------------------------------------------------------ */
  /* Dumpable                                                                 */
  /* ------------------------------------------------------------------------ */
public:
  void addDumpFieldToDumper(const std::string & dumper_name,
			    const std::string & field_id) override;

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
  
  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */

protected:
  Real time = 0;
  Real stepNum = 0;
  // first coefficient of friction
  SynchronizedArray<Real> mu_0;

  // coefficient of second term
  SynchronizedArray<Real> a;

  // coefficient of third term
  SynchronizedArray<Real> b;

  // velocity coefficient
  SynchronizedArray<Real> v_star;
  
  // state coefficient
  SynchronizedArray<Real> phi_star;

  // Dc the characteristic length
  SynchronizedArray<Real> d_c;

  // Revised R7S parameter
  SynchronizedArray<Real> v0;
  
};

/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

/// standard output stream operator
template <class Regularisation>
inline std::ostream & operator <<(std::ostream & stream, 
				  const NTNFricLawRateAndStateDieterich<Regularisation> & _this)
{
  _this.printself(stream);
  return stream;
}

} // namespace akantu

#include "ntn_friclaw_rate_and_state_dieterich_tmpl.hh"

#endif /* __AST_NTN_FRICLAW_RATE_AND_STATE_DIETERICH_HH__ */
