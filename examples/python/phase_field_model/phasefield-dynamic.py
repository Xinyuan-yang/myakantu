#!/usr/bin/env python
# coding: utf-8
__copyright__ = (
    "Copyright (©) 2021-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)"
    "Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)"
)
__license__ = "LGPLv3"


import akantu as aka

# reading material file
aka.parseInput("material.dat")
# creating mesh
spatial_dimension = 2
mesh = aka.Mesh(spatial_dimension)
mesh.read("plate.msh")


model = aka.CouplerSolidPhaseField(mesh)

solid = model.getSolidMechanicsModel()
phase = model.getPhaseFieldModel()

# initializing the Solid Mechanics Model with implicit solver for static
# resolution
solid.initFull(_analysis_method=aka._static)
solver = solid.getNonLinearSolver("static")
solver.set("max_iterations", 100)
solver.set("threshold", 1e-10)
solver.set("convergence_type", aka.SolveConvergenceCriteria.residual)

# adding another solver dynamic/quasi-static resolution (explicit Newmark with
# lumped mass)
solid.initNewSolver(aka._explicit_lumped_mass)

# initializing the PhaseField Model with linear implicit solver for static
# resolution
phase.initFull(_analysis_method=aka._static)

# initializing the PhaseField Model with Newton Raphson implicit solver for
# static resolution
phase.getNewSolver(
    "nonlinear_static",
    aka.TimeStepSolverType.static,
    aka.NonLinearSolverType.newton_raphson,
)
phase.setIntegrationScheme(
    "nonlinear_static", "damage", aka.IntegrationSchemeType.pseudo_time
)

solver = phase.getNonLinearSolver("nonlinear_static")
solver.set("max_iterations", 100)
solver.set("threshold", 1e-3)
solver.set("convergence_type", aka.SolveConvergenceCriteria.solution)


# Initialization for bulk vizualisation
solid.setBaseName("plate")
solid.addDumpFieldVector("displacement")
solid.addDumpFieldVector("external_force")
solid.addDumpFieldVector("velocity")
solid.addDumpField("strain")
solid.addDumpField("stress")
solid.addDumpField("damage")
solid.addDumpField("blocked_dofs")

# Dirichlet BC
solid.applyBC(aka.FixedValue(0.0, aka._x), "top")
solid.applyBC(aka.FixedValue(0.0, aka._x), "bottom")

solid.applyBC(aka.FixedValue(0.0, aka._x), "left")
solid.applyBC(aka.FixedValue(0.0, aka._x), "right")

# Pre-strain
solid.applyBC(aka.FixedValue(0.06e-3, aka._y), "top")
solid.applyBC(aka.FixedValue(-0.06e-3, aka._y), "bottom")

# Apply pre-strain by solving static problem
solid.solveStep("static")
solid.dump()


# #### **Damped dynamics resolution**
solid.setTimeStep(solid.getStableTimeStep() * 0.8)

# set maximum number of iteration
maxsteps = 1000
# solve using staggered scheme
for i in range(0, maxsteps):
    if i % 100 == 0:
        print("step {0}/{1}".format(i, maxsteps))
    model.solve("explicit_lumped", "")
    if i % 100 == 0:
        model.dump()
