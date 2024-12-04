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

namespace akantu{
  
  /* -------------------------------------------------------------------------- */
  struct StateEvolution {
    StateEvolution(){};
    virtual ~StateEvolution(){};
    virtual inline Real computeStateIncrement(Real slip_velocity, Real state, Real d_c, Real v0, Real time_step) = 0;
  };
  
  /* -------------------------------------------------------------------------- */
  struct AgingStateEvolution : StateEvolution {
    AgingStateEvolution(){};
    virtual ~AgingStateEvolution(){};
    virtual inline Real computeStateIncrement(Real slip_velocity, Real state, Real d_c, Real v0, Real time_step){
      return (1 - state*slip_velocity/d_c)*time_step;
    }
  };
  
/* -------------------------------------------------------------------------- */
struct SlipStateEvolution : StateEvolution {
  SlipStateEvolution(){};
  virtual ~SlipStateEvolution(){};
  virtual inline Real computeStateIncrement(Real slip_velocity, Real state, Real d_c, Real v0, Real time_step){
    return -(slip_velocity*state/d_c*log(slip_velocity*state/d_c))*time_step;
  }
};

/* -------------------------------------------------------------------------- */
struct RegularizedAgingStateEvolution : StateEvolution {
  RegularizedAgingStateEvolution(){};
  virtual ~RegularizedAgingStateEvolution(){};
  virtual inline Real computeStateIncrement(Real slip_velocity, Real state, Real d_c, Real v0, Real time_step){
    return (1-state*sqrt(slip_velocity*slip_velocity+v0*v0)/d_c)*time_step;
  }
};  


/* -------------------------------------------------------------------------- */
struct RandSFormulation {
  RandSFormulation(){};
  virtual ~RandSFormulation(){};
  virtual inline Real computeFrictionCoefficient(Real slip_velocity, Real state, Real mu_0, Real a, Real b, Real v0, Real v_star, Real state_star) = 0;
};
  
/* -------------------------------------------------------------------------- */
struct DieterichRandSFormulation : RandSFormulation{
  DieterichRandSFormulation(){};
  virtual ~DieterichRandSFormulation(){};
  virtual inline Real computeFrictionCoefficient(Real slip_velocity, Real state, Real mu_0, Real a, Real b, Real v0, Real v_star, Real state_star){
    return mu_0 + a * log(1+slip_velocity/v_star)+b*mu_0*log(1+state/state_star);
  }
};
  
} // namespace akantu
