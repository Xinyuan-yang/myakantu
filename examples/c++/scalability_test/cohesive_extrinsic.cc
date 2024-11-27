/**
 * @file   cohesive_extrinsic.cc
 *
 * @author Zineb Fouad <zineb.fouad@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author Seyedeh Mohadeseh Taheri Mousavi <mohadeseh.taherimousavi@epfl.ch>
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 *
 * @date creation: Tue May 08 2012
 * @date last modification: Wed Feb 06 2019
 *
 * @brief  Cohesive element examples in extrinsic
 *
 *
 * @section LICENSE
 *
 * Copyright (©) 2015-2021 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
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
 *
 */

/* -------------------------------------------------------------------------- */
#include <non_linear_solver.hh>
#include <solid_mechanics_model_cohesive.hh>
/* -------------------------------------------------------------------------- */
#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
/* -------------------------------------------------------------------------- */

using clk = std::chrono::high_resolution_clock;
using seconds = std::chrono::duration<double>;
using milliseconds = std::chrono::duration<double, std::milli>;

// #define AKANTU_VERSION_MAJOR 2
class Chrono {
public:
  Chrono(int prank, int psize) : prank(prank), psize(psize) {}

  inline void start() { _start = clk::now(); };
  inline void store_time(const std::string & type) {
    clk::time_point _end = clk::now();
    if (measures.find(type) == measures.end()) {
      measures[type] = _end - _start;
      nb_measures[type] = 1;
      if (prank == 0) {
        std::cout << "Passing the first " << type << " chrono! ["
                  << measures[type].count() << "]" << std::endl;
      }
    } else {
      measures[type] += _end - _start;
      ++nb_measures[type];
    }

    _start = clk::now();
  }

  virtual void printself(std::ostream & stream, int indent = 0) const {
    std::string space(indent, AKANTU_INDENT);

    stream << space << "Chrono [" << std::endl;
    for (auto && measure : measures) {
      const unsigned int & nb_measure = nb_measures.find(measure.first)->second;
      stream << space << " + " << measure.first << "\t: " << std::setw(25)
             << std::fixed << std::setprecision(16) << measure.second.count()
             << "us - nb_repetition: " << nb_measure << std::endl;
    }
    stream << space << "]" << std::endl;
  }

  virtual void printself_csv(std::ostream & stream, int indent = 0) const {
    std::string space(indent, AKANTU_INDENT);
    stream << "\"psize\"";
    for (auto && measure : measures) {
      stream << ", \"" << measure.first << "\""
             << ", \"" << measure.first << " nb_rep"
             << "\"";
    }
    stream << std::endl;

    stream << psize;
    for (auto && measure : measures) {
      const unsigned int & nb_measure = nb_measures.find(measure.first)->second;
      stream << ", " << measure.second.count() << ", " << nb_measure;
    }
    stream << std::endl;
  }

private:
  clk::time_point _start;
  std::map<std::string, seconds> measures;
  std::map<std::string, unsigned int> nb_measures;
  int prank, psize;
};

inline std::ostream & operator<<(std::ostream & stream, const Chrono & _this) {
  _this.printself_csv(stream);
  return stream;
}

using namespace akantu;

int main(int argc, char * argv[]) {
  initialize("material-elastic.dat", argc, argv);

  const UInt spatial_dimension = 3;

  auto prank = Communicator::getStaticCommunicator().whoAmI();
  auto psize = Communicator::getStaticCommunicator().getNbProc();
  Chrono chrono(prank, psize);

  const auto & usersect = getUserParser();

  const Real c = usersect.getParameter("compression");
  const Real s = usersect.getParameter("shear");
  const Real inc_s = usersect.getParameter("inc_shear");
  const bool output_energy = usersect.getParameter("output_energy", true);
  const bool output_paraview = usersect.getParameter("output_paraview", true);
  const bool cohesive_insertion =
      usersect.getParameter("cohesive_insertion", true);
  const UInt max_steps = usersect.getParameter("max_steps");
  const std::string mesh_filename = usersect.getParameter("mesh");

  if (prank == 0) {
    std::cout << "Paramters:\n"
              << " - output_energy: " << output_energy << "\n"
              << " - output_paraview: " << output_paraview << "\n"
              << " - cohesive_insertion: " << cohesive_insertion << "\n"
              << " - max_steps: " << max_steps << "\n"
              << " - mesh_filename: " << mesh_filename << "\n";
  }
  chrono.start();

  clk::time_point start_time = clk::now();

  Mesh mesh(spatial_dimension);

  if (prank == 0) {
    mesh.read(mesh_filename);
    chrono.store_time("read_mesh");
  }
  mesh.distribute();

  chrono.store_time("dist_mesh");

  SolidMechanicsModelCohesive model(mesh);

  /// model initialization
  model.initFull(_analysis_method = _static, _is_extrinsic = true);

  chrono.store_time("init_full");

  auto & blocked_dofs = model.getBlockedDOFs();
  auto & force = model.getExternalForce();

  /// boundary conditions
  model.applyBC(BC::Dirichlet::FixedValue(0.0, _z), "bottom");     // face
  model.applyBC(BC::Dirichlet::FixedValue(0.0, _x), "right line"); // line
  blocked_dofs(3, _y) = true;                                      // point

  Matrix<Real> compression{{c, 0., 0.}, {0., c, 0.}, {0., 0., c}};
  Matrix<Real> shear{{0., 0., s}, {0., 0., 0.}, {s, 0., 0.}};

  force.zero();

  model.applyBC(BC::Neumann::FromHigherDim(compression), "top");
  model.applyBC(BC::Neumann::FromHigherDim(compression), "bottom");

  model.applyBC(BC::Neumann::FromHigherDim(shear), "top");
  model.applyBC(BC::Neumann::FromHigherDim(shear), "bottom");
  model.applyBC(BC::Neumann::FromHigherDim(shear), "side left");
  model.applyBC(BC::Neumann::FromHigherDim(shear), "side right");

  if (output_paraview) {
    model.setBaseName("extrinsic");
    model.addDumpFieldVector("displacement");
    model.addDumpField("internal_force");
    model.addDumpField("external_force");
    model.addDumpField("stress");
    model.addDumpField("blocked_dofs");
    model.addDumpField("grad_u");

    model.addDumpFieldToDumper("cohesive elements", "displacement");
    model.addDumpFieldToDumper("cohesive elements", "tractions");

    model.dump();
    model.dump("cohesive elements");
  }

  chrono.store_time("initial_conditons");

  model.solveStep("static");
  chrono.store_time("static_solve");

  model.initNewSolver(_explicit_lumped_mass);

  std::ofstream fout;
  if (output_energy and prank == 0) {
    fout.open("energies.csv", std::ofstream::out | std::ofstream::trunc);
    fout << "step, ed, ep, ek, ew, et" << std::endl;
  }

  Real Ed{0}, Ep{0}, Ek{0}, Ew{0};
  if (output_energy) {
    Ep = model.getEnergy("potential");
    Ek = model.getEnergy("kinetic");
    Ew += model.getEnergy("external work");
  }
  auto Et = Ed + Ep + Ek - Ew;

  if (output_energy and prank == 0) {
    fout << 0 << ", " << Ed << ", " << Ep << ", " << Ek << ", " << Ew << ", "
         << Et << std::endl;
  }

  Real time_step = model.getStableTimeStep() * 0.05;
  model.setTimeStep(time_step);
  if (prank == 0) {
    std::cout << "Time step: " << time_step << std::endl;
  }

  if (output_paraview) {
    model.addDumpField("velocity");
    model.addDumpField("acceleration");
    model.addDumpFieldToDumper("cohesive elements", "velocity");

    model.dump();
    model.dump("cohesive elements");
  }

  Matrix<Real> new_shear{{0., 0., inc_s}, {0., 0., 0.}, {inc_s, 0., 0.}};

  seconds init_time = clk::now() - start_time;
  chrono.store_time("before_step");

  auto nb_dofs_start = mesh.getNbGlobalNodes() * spatial_dimension;

  start_time = clk::now();
  /// Main loop
  for (auto s : arange(1, max_steps + 1)) {
    if (s % 100 == 0 and s < 10000) {
      model.applyBC(BC::Neumann::FromHigherDim(new_shear), "top");
      model.applyBC(BC::Neumann::FromHigherDim(new_shear), "bottom");
      model.applyBC(BC::Neumann::FromHigherDim(new_shear), "side left");
      model.applyBC(BC::Neumann::FromHigherDim(new_shear), "side right");
      chrono.store_time("boundary_conditions");
    }
    if (cohesive_insertion) {
      model.checkCohesiveStress();
      chrono.store_time("check_cohesive_stress");
    }
    model.solveStep("explicit_lumped");
    chrono.store_time("solve_step");

    if (output_energy) {
      Ed = model.getEnergy("dissipated");
      Ep = model.getEnergy("potential");
      Ek = model.getEnergy("kinetic");
      Ew += model.getEnergy("external work");

      Et = Ed + Ep + Ek - Ew;

      if (prank == 0) {
        fout << s << ", " << Ed << ", " << Ep << ", " << Ek << ", " << Ew
             << ", " << Et << std::endl;
      }
      chrono.store_time("energies");
    }

    if (output_paraview and (s % 100 == 0)) {
      model.dump();
      model.dump("cohesive elements");

      if (prank == 0) {
        milliseconds loop_time = clk::now() - start_time;
        std::cout << "passing step " << s << "/" << max_steps
                  << " - nb_cohesive_element: "
                  << mesh.getNbElement(spatial_dimension, _not_ghost,
                                       _ek_cohesive)
                  << " - nb_dofs: "
                  << mesh.getNbGlobalNodes() * spatial_dimension << " - "
                  << loop_time.count() / s << "\t\t \r";
        std::cout.flush();
      }
      chrono.store_time("dumpers");
    }
  }
  auto nb_dofs_end = mesh.getNbGlobalNodes() * spatial_dimension;

  seconds loop_time = clk::now() - start_time;
  auto nb_cohesive_elements =
      mesh.getNbElement(spatial_dimension, _not_ghost, _ek_cohesive);
  Communicator::getStaticCommunicator().allReduce(nb_cohesive_elements);
  Ed = model.getEnergy("dissipated");
  if (prank == 0) {
    std::cout << std::endl;
    std::cout << "Cohesive info: dissipated energy: " << Ed
              << " - nb_cohesive_element: " << nb_cohesive_elements
              << std::endl;
    std::cout << "Nb proc: "
              << Communicator::getStaticCommunicator().getNbProc() << std::endl;
    std::cout << "Full time: " << (init_time + loop_time).count() << std::endl;
    std::cout << "Init time: " << init_time.count() << std::endl;
    std::cout << "Step time: " << loop_time.count()
              << " - nb_steps: " << max_steps
              << " - time_per_step: " << (loop_time.count() / max_steps)
              << std::endl;
    std::cout << "Ns DOFs - start: " << nb_dofs_start
              << " - end: " << nb_dofs_end << std::endl;
  }

  if (prank == 0) {
    std::cerr << chrono << std::endl;
  }

  if (output_energy and prank == 0) {
    fout.close();
  }

  finalize();

  return EXIT_SUCCESS;
}
