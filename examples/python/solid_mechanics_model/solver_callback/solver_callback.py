#!/usr/bin/env python3
__copyright__ = (
    "Copyright (©) 2021-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)"
    "Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)"
)
__license__ = "LGPLv3"


import numpy as np
import akantu as aka


class SolverCallback(aka.InterceptSolverCallback):
    def __init__(self, model):
        super().__init__(model)
        self.model = model

        mesh = model.getMesh()
        left = mesh.getElementGroup("Left").getNodeGroup().getNodes()
        right = mesh.getElementGroup("Right").getNodeGroup().getNodes()
        position = mesh.getNodes()

        self.pair = []
        for node_l in left:
            node_l = int(node_l)
            for node_r in right:
                node_r = int(node_r)
                if abs(position[node_r, 1] - position[node_l, 1]) < 1e-6:
                    self.pair.append([node_l, node_r])

        self.periodic_K_modif = aka.TermsToAssemble(
            "displacement", "displacement"
        )

        matrix_type = self.model.getMatrixType("K")
        for p in self.pair:
            # blocked_dofs[p[1]] = True
            # a u_{i, x} + b u_{j, x} = 0
            # self.periodic_K_modif(i*dim + aka._x, i*dim + aka._x, a)
            # self.periodic_K_modif(i*dim + aka._x, j*dim + aka._x, b)

            self.periodic_K_modif(p[0] * 2, p[0] * 2, 1)
            self.periodic_K_modif(p[0] * 2, p[1] * 2, -1)
            if matrix_type == aka._unsymmetric:
                self.periodic_K_modif(p[1] * 2, p[0] * 2, -1)
                self.periodic_K_modif(p[1] * 2, p[1] * 2, 1)

        self.first = True
        self.k_release = -1

    def assembleMatrix(self, matrix_id):
        self.model.assembleMatrix(matrix_id)
        if matrix_id == "K":
            release = self.model.getDOFManager().getMatrix("K").getRelease()
            if release == self.k_release:
                return

            if self.first:
                self.model.getDOFManager().getMatrix("K").saveMatrix("K0.mtx")

            self.model.getDOFManager().assemblePreassembledMatrix(
                "K", self.periodic_K_modif
            )

            if self.first:
                self.model.getDOFManager().getMatrix("K").saveMatrix("K1.mtx")

            self.k_release = (
                self.model.getDOFManager().getMatrix("K").getRelease()
            )
            self.first = False

    def assembleResidual(self):
        displacement = self.model.getDisplacement()
        force = np.zeros(displacement.shape)
        for p in self.pair:
            force[p[0], 0] += displacement[p[0], 0] - displacement[p[1], 0]
            force[p[1], 0] += displacement[p[1], 0] - displacement[p[0], 0]
        self.model.getDOFManager().assembleToResidual(
            "displacement", force, -1.0
        )
        self.model.assembleResidual()


# -----------------------------------------------------------------------------
def main():
    spatial_dimension = 2
    mesh_file = "bar.msh"
    max_steps = 250
    time_step = 1e-3

    aka.parseInput("material.dat")

    # -------------------------------------------------------------------------
    # Initialization
    # -------------------------------------------------------------------------
    mesh = aka.Mesh(spatial_dimension)
    mesh.read(mesh_file)

    model = aka.SolidMechanicsModel(mesh)

    model.initFull(_analysis_method=aka._implicit_dynamic)

    model.setBaseName("solver_callback")
    model.addDumpFieldVector("displacement")
    model.addDumpFieldVector("acceleration")
    model.addDumpFieldVector("velocity")
    model.addDumpFieldVector("internal_force")
    model.addDumpFieldVector("external_force")
    model.addDumpField("strain")
    model.addDumpField("stress")
    model.addDumpField("blocked_dofs")

    # -------------------------------------------------------------------------
    # boundary conditions
    # -------------------------------------------------------------------------
    model.applyBC(aka.FixedValue(0, aka._y), "YBlocked")

    # -------------------------------------------------------------------------
    # initial conditions
    # -------------------------------------------------------------------------
    displacement = model.getDisplacement()
    velocity = model.getVelocity()
    position = mesh.getNodes()

    L = 1  # pulse_width
    A = 0.01
    v = np.sqrt(
        model.getMaterial(0).getReal("E") / model.getMaterial(0).getReal("rho")
    )
    k = 0.1 * 2 * np.pi * 3 / L
    t = 0.0
    velocity[:, 0] = k * v * A * np.sin(k * ((position[:, 0] - 5.0) - v * t))
    displacement[:, 0] = A * np.cos(k * ((position[:, 0] - 5.0) - v * t))

    # -------------------------------------------------------------------------
    # timestep value computation
    # -------------------------------------------------------------------------
    time_factor = 0.8
    stable_time_step = model.getStableTimeStep() * time_factor

    print("Stable Time Step = {0}".format(stable_time_step))
    print("Required Time Step = {0}".format(time_step))

    time_step = stable_time_step * time_factor

    model.setTimeStep(time_step)
    solver_callback = SolverCallback(model)

    solver = model.getNonLinearSolver()
    solver.set("max_iterations", 100)
    solver.set("threshold", 1e-7)

    # -------------------------------------------------------------------------
    # loop for evolution of motion dynamics
    # -------------------------------------------------------------------------
    print("step,step * time_step,epot,ekin,epot + ekin")
    for step in range(0, max_steps + 1):
        model.solveStep(solver_callback)
        # model.solveStep()

        if step % 10 == 0:
            model.dump()

        epot = model.getEnergy("potential")
        ekin = model.getEnergy("kinetic")

        # output energy calculation to screen
        print(
            "{0},{1},{2},{3},{4}".format(
                step, step * time_step, epot, ekin, (epot + ekin)
            )
        )

    return


# -----------------------------------------------------------------------------
if __name__ == "__main__":
    main()
