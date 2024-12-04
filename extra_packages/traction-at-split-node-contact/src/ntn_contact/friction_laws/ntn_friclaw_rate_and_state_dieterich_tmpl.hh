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
#include "dumper_text.hh"

namespace akantu{
  
/* -------------------------------------------------------------------------- */
template <class Regularisation>
NTNFricLawRateAndStateDieterich<Regularisation>::NTNFricLawRateAndStateDieterich(NTNBaseContact & contact,
										 const ID & id) :
  NTNFricLawCoulomb<Regularisation>(contact,id),
  mu_0(0,1,0.,id+":mu_0",0.,"mu_0"),
  a(0,1,0.,id+":a",0.,"a"),
  b(0,1,0.,id+":b",0.,"b"),
  v_star(0,1,0.,id+":v_star",0.,"v_star"),
  phi_star(0,1,0.,id+":phi_star",0.,"phi_star"),
  d_c(0,1,0.,id+":d_c",0.,"d_c"),
  v0(0,1,0.,id+":v0",0.,"v0"){
  AKANTU_DEBUG_IN();
  
  NTNFricLawCoulomb<Regularisation>::registerSynchronizedArray(this->mu_0);
  NTNFricLawCoulomb<Regularisation>::registerSynchronizedArray(this->a);
  NTNFricLawCoulomb<Regularisation>::registerSynchronizedArray(this->b);
  NTNFricLawCoulomb<Regularisation>::registerSynchronizedArray(this->v_star);
  NTNFricLawCoulomb<Regularisation>::registerSynchronizedArray(this->phi_star);
  NTNFricLawCoulomb<Regularisation>::registerSynchronizedArray(this->d_c);
  NTNFricLawCoulomb<Regularisation>::registerSynchronizedArray(this->v0);
  
  this->registerParam("mu_0",           this->mu_0, _pat_parsmod, "first coefficient of friction");
  this->registerParam("a",              this->a, _pat_parsmod, "coefficient of second term");
  this->registerParam("b",              this->b,  _pat_parsmod, "coefficient of third term");
  this->registerParam("v_star",         this->v_star, _pat_parsmod, "velocity coefficient");
  this->registerParam("phi_star",       this->phi_star, _pat_parsmod, "state coefficient");
  this->registerParam("friction_state", this->friction_state, _pat_parsmod, "state of friction");
  this->registerParam("d_c",            this->d_c,  _pat_parsmod, "characteristic length");
  this->registerParam("v0",             this->v0,  _pat_parsmod, "R&S parameter");
  
  this->setParameterAccessType("mu", _pat_readable);
  
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <class Regularisation>
void NTNFricLawRateAndStateDieterich<Regularisation>::computeFrictionalStrength() {
  AKANTU_DEBUG_IN();
  
  computeFrictionCoefficient();
  NTNFricLawCoulomb<Regularisation>::computeFrictionalStrength();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <class Regularisation>
void NTNFricLawRateAndStateDieterich<Regularisation>::computeFrictionCoefficient() {
  AKANTU_DEBUG_IN();
  
  // get arrays
  const SynchronizedArray<bool> & stick = this->internalGetIsSticking();
  const SynchronizedArray<Real> & slip_velocity_norm  = this->internalGetSlipVelocity();
  
  // Retrieve time step to calculate current time
  Real time_step = this->contact.getModel().getTimeStep();
  time += time_step;
  stepNum += 1; 

  UInt nb_contact_nodes = this->contact.getNbContactNodes();
  for (UInt n=0; n<nb_contact_nodes; ++n) {
    //=========================================================================================================================================================================     
    // Dieterich rate and state law
    // With regularized evolution law !
    // A + 1 term is added inside the log of velocity to avoid NaN values for v=0
    // sticking phase
    if (stick(n)) {
      this->mu(n) = this->mu_0(n) + this->a(n) * log(1+slip_velocity_norm(n)/this->v_star(n)) + this->b(n) * this->mu_0(n) * log(this->friction_state(n)/this->phi_star(n));
    }
    // slipping phase
    else {
      this->friction_state(n) += (1 - this->friction_state(n)*sqrt(pow(slip_velocity_norm(n),2)+pow(this->v0(n),2))/this->d_c(n))*time_step;
      this->mu(n) = this->mu_0(n) + this->a(n) * log(1+slip_velocity_norm(n)/this->v_star(n)) + this->b(n) * this->mu_0(n) * log(this->friction_state(n)/this->phi_star(n));
    }
    //=========================================================================================================================================================================    
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <class Regularisation>
void NTNFricLawRateAndStateDieterich<Regularisation>::addNoiseToStateVariable(std::map<Real,Real> &noise) {
  AKANTU_DEBUG_IN();
  
  std::cout<< "Adding noise to the state field" << std::endl;
  
  const Array<Real> & cur_pos = this->contact->getModel().getCurrentPosition();
  const SynchronizedArray<UInt> & slaves  = this->internalGetSlaves();
  UInt nb_contact_nodes = this->contact->getNbContactNodes();
  
  for (UInt n=0; n<nb_contact_nodes; ++n) {
    UInt slave = slaves(n);

    auto it = noise.lower_bound(cur_pos(slave));
    
    Real perturbation = it->second;
    
    std::cout << "Slave :" << slave << " " << cur_pos(slave) << " : perturbation is " << perturbation << std::endl;

    this->friction_state(n) = this->friction_state(n) * (1 + perturbation);
  }
   
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <class Regularisation>
void NTNFricLawRateAndStateDieterich<Regularisation>::registerSynchronizedArray(SynchronizedArrayBase & array) {
  AKANTU_DEBUG_IN();

  this->mu_0.registerDependingArray(array);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <class Regularisation>
void NTNFricLawRateAndStateDieterich<Regularisation>::dumpRestart(const std::string & file_name) const {
  AKANTU_DEBUG_IN();

  this->mu_0.dumpRestartFile(file_name);
  this->a.dumpRestartFile(file_name);
  this->b.dumpRestartFile(file_name);
  this->phi_star.dumpRestartFile(file_name);
  this->v_star.dumpRestartFile(file_name);
  this->d_c.dumpRestartFile(file_name);
  this->v0.dumpRestartFile(file_name);
  this->friction_state.dumpRestartFile(file_name);

  NTNFricLawCoulomb<Regularisation>::dumpRestart(file_name);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <class Regularisation>
void NTNFricLawRateAndStateDieterich<Regularisation>::readRestart(const std::string & file_name) {
  AKANTU_DEBUG_IN();

  this->mu_0.readRestartFile(file_name);
  this->a.readRestartFile(file_name);
  this->b.readRestartFile(file_name);
  this->phi_star.readRestartFile(file_name);
  this->v_star.readRestartFile(file_name);
  this->d_c.readRestartFile(file_name);
  this->v0.readRestartFile(file_name);
  this->friction_state.readRestartFile(file_name);

  NTNFricLawCoulomb<Regularisation>::readRestart(file_name);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <class Regularisation>
void NTNFricLawRateAndStateDieterich<Regularisation>::printself(std::ostream & stream, int indent) const {
  AKANTU_DEBUG_IN();
  std::string space;
  for(Int i = 0; i < indent; i++, space += AKANTU_INDENT);

  stream << space << "NTNFricLawRateAndStateDieterich [" << std::endl;
  NTNFricLawCoulomb<Regularisation>::printself(stream, ++indent);
  stream << space << "]" << std::endl;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <class Regularisation>
void NTNFricLawRateAndStateDieterich<Regularisation>::addDumpFieldToDumper(const std::string & dumper_name,
								  const std::string & field_id) {
  AKANTU_DEBUG_IN();

  if(field_id == "mu_0") {
    this->internalAddDumpFieldToDumper(dumper_name,
				       field_id,
				       std::make_unique<dumpers::NodalField<Real>>(this->mu_0.getArray()));
  }
  else if(field_id == "a") {
    this->internalAddDumpFieldToDumper(dumper_name,
				       field_id,
				       std::make_unique<dumpers::NodalField<Real>>(this->a.getArray()));
  }
  else if(field_id == "b") {
    this->internalAddDumpFieldToDumper(dumper_name,
				       field_id,
				       std::make_unique<dumpers::NodalField<Real>>(this->b.getArray()));
  }
  else if(field_id == "phi_star") {
    this->internalAddDumpFieldToDumper(dumper_name,
				       field_id,
				       std::make_unique<dumpers::NodalField<Real>>(this->phi_star.getArray()));
  }
  else if(field_id == "v_star") {
    this->internalAddDumpFieldToDumper(dumper_name,
				       field_id,
				       std::make_unique<dumpers::NodalField<Real>>(this->v_star.getArray()));
  }
  else if(field_id == "friction_state") {
    this->internalAddDumpFieldToDumper(dumper_name,
				       field_id,
				       std::make_unique<dumpers::NodalField<Real>>(this->friction_state.getArray()));
  }   
  else if(field_id == "d_c") {
    this->internalAddDumpFieldToDumper(dumper_name,
				       field_id,
				       std::make_unique<dumpers::NodalField<Real>>(this->d_c.getArray()));
  }
  else {
    NTNFricLawCoulomb<Regularisation>::addDumpFieldToDumper(dumper_name, field_id);
  }

  AKANTU_DEBUG_OUT();
}

} // namespace akantu
