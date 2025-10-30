#!/usr/bin/env python3
__copyright__ = (
    "Copyright (©) 2018-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)"
    "Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)"
)
__license__ = "LGPLv3"


import sys
from patch_test_linear_heat_transfer_fixture import TestPatchTestHTMLinear
import akantu


def foo(self):
    self.initModel(akantu._static, "heat_transfer_input.dat")

    solver = self.model.getNonLinearSolver()
    solver.set("max_iterations", 2)
    solver.set("threshold", 2e-4)
    solver.set("convergence_type", akantu.SolveConvergenceCriteria.residual)

    self.model.solveStep()
    self.checkAll()


def test():
    TestPatchTestHTMLinear.TYPED_TEST(foo, "Static")


if "pytest" not in sys.modules:
    test()
