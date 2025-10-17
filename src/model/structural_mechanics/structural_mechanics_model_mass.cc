/**
 * Copyright (©) 2014-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "aka_math.hh"
#include "integrator_gauss.hh"
#include "material.hh"
#include "shape_structural.hh"
#include "structural_mechanics_model.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {

class ComputeRhoFunctorStruct {
public:
  explicit ComputeRhoFunctorStruct(const StructuralMechanicsModel & model)
      : model(model) {};

  void operator()(Matrix<Real> & rho, const Element & element) const {
    const auto & mat = model.getMaterialByElement(element);
    const Real mat_rho = mat.rho; // This is the density
    const Real mat_A = mat.A;     // This is the cross section
    const Real mat_mass_per_length =
        mat_rho * mat_A; // Mass of the beam per unit length

    rho.set(
        mat_mass_per_length); // The integrator _recquiers_ mass per unit length
  }

private:
  const StructuralMechanicsModel & model;
};

/* -------------------------------------------------------------------------- */
void StructuralMechanicsModel::assembleMassMatrix() {
  AKANTU_DEBUG_IN();

  if (not need_to_reassemble_mass) {
    return;
  }

  if (not getDOFManager().hasMatrix("M")) {
    getDOFManager().getNewMatrix("M", getMatrixType("M"));
  }

  this->getDOFManager().zeroMatrix("M");
  assembleMassMatrix(_not_ghost);

  need_to_reassemble_mass = false;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void StructuralMechanicsModel::assembleMassMatrix(GhostType ghost_type) {
  AKANTU_DEBUG_IN();
  auto & fem = getFEEngineClass<MyFEEngineType>();
  ComputeRhoFunctorStruct compute_rho(*this);

  for (auto type :
       mesh.elementTypes(spatial_dimension, ghost_type, _ek_structural)) {
    fem.assembleFieldMatrix(compute_rho, "M", "displacement",
                            this->getDOFManager(), type, ghost_type);
  }
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void StructuralMechanicsModel::assembleLumpedMassMatrix() {

  if (not this->need_to_reassemble_lumped_mass) {
    return;
  }

  allocNodalField(this->mass, nb_degree_of_freedom, "lumped_mass");

  if (!this->getDOFManager().hasLumpedMatrix("M")) {
    this->getDOFManager().getNewLumpedMatrix("M");
  }

  this->getDOFManager().zeroLumpedMatrix("M");

  // Preparing
  const UInt nb_nodes = mesh.getNbNodes();
  const auto & nodes = mesh.getNodes();
  const auto & node_it = make_view(nodes, spatial_dimension).begin();
  auto & lumped_mass = *(this->mass);

  lumped_mass.zero(); // Set the lumped mass to zero

  Array<Real> volumes(nb_nodes, 1, 0., "volumes");

  /// We now compute the mass and the volume, but not the inertia
  for (auto type :
       mesh.elementTypes(spatial_dimension, _not_ghost, _ek_structural)) {
    const auto & element_material_id = this->element_material(type);
    const auto & connectivity = mesh.getConnectivity(type);
    // const auto nb_element = connectivity.size();

    if (type != _bernoulli_beam_3 and type != _bernoulli_beam_2) {
      AKANTU_EXCEPTION(
          "The lumped mass was not implemented for non Bernoulli Beams");
    }

    // Now iterate through all the connectivity
    for (auto && data : zip(make_view(connectivity, 2), element_material_id)) {
      const auto & conn = std::get<0>(data);
      auto node1 = conn(0);
      auto node2 = conn(0);

      const auto & material = this->materials.at(std::get<1>(data));
      const auto rho = material.rho;
      const auto cross_section = material.A;

      auto && coord_node_1 = node_it[node1];
      auto && coord_node_2 = node_it[node2];

      const auto length = coord_node_1.distance(coord_node_2);
      const auto volume = length * cross_section;
      const auto mass = volume * rho;

      // Now distribute the mass at the right places
      for (auto d : arange(spatial_dimension)) {
        lumped_mass(node1, d) += mass / 2.;
        lumped_mass(node2, d) += mass / 2.;
      }; // end for(d):

      volumes(node1) += volume / 2.; // this entry never point to mass
      volumes(node2) += volume / 2.;
    }
  }

  const auto pi = std::atan(1.) * 4;

  // We now compute the inertia.
  if (spatial_dimension == 2) {
    /* This is the 2D case, so we are assuming that the mass is inside a disc.
     * Which is given as
     * \begin{align}
     *   I := m \cdot r^2
     * \end{align}
     *
     * The radius is obtained by assuming that the volume, that is associated to
     * a beam, forms a uniformly disc. From this volume we can compute the
     * radius.
     */
    for (auto && data :
         zip(make_view(lumped_mass, nb_degree_of_freedom), volumes)) {
      const Real volume = std::get<1>(data);
      auto & masses = std::get<0>(data);

      const Real r2 = volume / pi; // The square of the volume

      const Real inertia = masses(0) * r2;
      masses(spatial_dimension) = inertia;
    }
  } else if (spatial_dimension == 3) {
    /* This is essentially the same as for 2D only that we assume here,
     * that the mass is uniformly distributed inside a sphere.
     * And thus we have
     * \begin{align}
     *   I := \frac{2}{5} m \cdot r^2
     * \end{align}
     *
     * We also have to set three values, for the three axis.
     * But since we assume a sphere, they are all the same.
     */
    for (auto && data :
         zip(make_view(lumped_mass, nb_degree_of_freedom), volumes)) {
      const Real volume = std::get<1>(data);
      auto & masses = std::get<0>(data);

      const Real r2 = std::pow((volume * 3.) / (4 * pi), 2. / 3.);
      const Real inertia = (2. / 5.) * masses(0) * r2;

      for (auto && m : masses) {
        m = inertia;
      }
    }
  } else {
    AKANTU_EXCEPTION("The dimension " << spatial_dimension << " is not known.");
  }

  this->getDOFManager().assembleToLumpedMatrix("displacement", lumped_mass,
                                               "M");

  this->need_to_reassemble_lumped_mass = false;
}

} // namespace akantu
