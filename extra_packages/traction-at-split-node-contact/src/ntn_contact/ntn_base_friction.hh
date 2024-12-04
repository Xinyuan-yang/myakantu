/**
 * Copyright (©) 2016-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#ifndef AST_NTN_BASE_FRICTION_HH_
#define AST_NTN_BASE_FRICTION_HH_

/* -------------------------------------------------------------------------- */
// akantu
#include "parsable.hh"

// simtools
#include "ntn_base_contact.hh"

/* -------------------------------------------------------------------------- */
namespace akantu {

/* -------------------------------------------------------------------------- */
template <>
inline void ParameterTyped<akantu::SynchronizedArray<Real>>::setAuto(
    const ParserParameter & in_param) {
  Parameter::setAuto(in_param);
  Real r = in_param;
  param.setAndChangeDefault(r);
}

/* -------------------------------------------------------------------------- */
template <>
template <>
inline void ParameterTyped<akantu::SynchronizedArray<Real>>::setTyped<Real>(
    const Real & value) {
  param.setAndChangeDefault(value);
}

/* -------------------------------------------------------------------------- */
class NTNBaseFriction : public Parsable, public Dumpable {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  NTNBaseFriction(NTNBaseContact & contact, const ID & id = "friction");
  ~NTNBaseFriction() override = default;

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// compute friction traction
  virtual void computeFrictionTraction();

  /// assemble the friction traction
  virtual void assembleGlobalFrictionTraction();

  /// compute stick traction (friction traction needed to stick the nodes)
  virtual void computeStickTraction();

  /// compute slip
  virtual void updateSlip();

  /// register Syncronizedarrays for sync
  virtual void registerSynchronizedArray(SynchronizedArrayBase & array);

  /// dump restart file
  virtual void dumpRestart(const std::string & file_name) const;

  /// read restart file
  virtual void readRestart(const std::string & file_name);

  /// set to steady state
  virtual void setToSteadyState() { AKANTU_TO_IMPLEMENT(); };

  /// get the number of sticking nodes (in parallel)
  /// a node that is not in contact does not count as sticking
  virtual UInt getNbStickingNodes() const;

  /// function to print the contain of the class
  void printself(std::ostream & stream, int indent = 0) const override;

protected:
  /// compute frictional strength according to friction law
  virtual void computeFrictionalStrength() = 0;

  /* ------------------------------------------------------------------------ */
  /* Dumpable                                                                 */
  /* ------------------------------------------------------------------------ */
public:
  void addDumpFieldToDumper(const std::string & dumper_name,
                            const std::string & field_id) override;

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  AKANTU_GET_MACRO(Contact, contact, const NTNBaseContact &)

  AKANTU_GET_MACRO(IsSticking, is_sticking, const SynchronizedArray<bool> &)
  AKANTU_GET_MACRO(FrictionalStrength, frictional_strength,
                   const SynchronizedArray<Real> &)
  AKANTU_GET_MACRO(FrictionTraction, friction_traction,
                   const SynchronizedArray<Real> &)
  /// get the NTNBaseFriction::friction_traction array (friction_traction)
  AKANTU_GET_MACRO_DEREF_PTR_NOT_CONST(GlobalFrictionTraction,
                                       global_friction_traction)
  /// get the  NTNBaseFriction::friction_traction array (friction_traction)
  AKANTU_GET_MACRO_DEREF_PTR(GlobalFrictionTraction, global_friction_traction)
  AKANTU_GET_MACRO(Slip, slip, const SynchronizedArray<Real> &)
  AKANTU_GET_MACRO(CumulativeSlip, cumulative_slip,
                   const SynchronizedArray<Real> &)
  AKANTU_GET_MACRO(SlipVelocity, slip_velocity, const SynchronizedArray<Real> &)
  AKANTU_GET_MACRO(SlipVelocityNorm,    slip_velocity_norm, const SynchronizedArray<Real> &)
  //AKANTU_GET_MACRO(State,                              friction_state, const SynchronizedArray<Real> &)  

  /// set parameter of a given node
  /// (if you need to set to all: used the setMixed function of the Parsable).
  virtual void setParam(const std::string & name, UInt node, Real value);

  // replaced by the setMixed of the Parsable
  // virtual void setParam(const std::string & param, Real value) {
  //   AKANTU_ERROR("Friction does not know the following parameter: " <<
  //   param);
  // };

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  NTNBaseContact & contact;

  // if node is sticking
  SynchronizedArray<bool> is_sticking;
  // frictional strength
  SynchronizedArray<Real> frictional_strength;
  // friction force
  SynchronizedArray<Real> friction_traction;
  // slip
  SynchronizedArray<Real> slip;
  SynchronizedArray<Real> cumulative_slip;
  // slip velocity (tangential vector)
  SynchronizedArray<Real> slip_velocity;
  SynchronizedArray<Real> slip_velocity_norm;
  
  // friction force
  std::unique_ptr<Array<Real>> global_friction_traction;
};

/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

// #include "ntn_base_friction_inline_impl.hh"

/// standard output stream operator
inline std::ostream & operator<<(std::ostream & stream,
                                 const NTNBaseFriction & _this) {
  _this.printself(stream);
  return stream;
}

} // namespace akantu

#endif /* AST_NTN_BASE_FRICTION_HH_ */
