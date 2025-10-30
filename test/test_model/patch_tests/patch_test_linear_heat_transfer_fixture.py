#!/usr/bin/env python3
__copyright__ = (
    "Copyright (©) 2018-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)"
    "Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)"
)
__license__ = "LGPLv3"


import patch_test_linear_fixture
import akantu


class TestPatchTestHTMLinear(patch_test_linear_fixture.TestPatchTestLinear):
    model_type = akantu.HeatTransferModel

    def applyBC(self):
        super().applyBC()
        temperature = self.model.getTemperature()
        self.applyBConDOFs(temperature)

    def checkAll(self):
        temperature = self.model.getTemperature()
        cl = self.model.getConstitutiveLaw(0)
        C = cl.getMatrix("conductivity")

        self.checkDOFs(temperature)

        grad_t = cl.getInternalReal("∇u")
        self.checkGradient(grad_t(self.elem_type), temperature)

        self.prescribed_gradient(temperature)
        flux = cl.getInternalReal("D∇u")
        self.checkResults(
            lambda grad_T: C.dot(grad_T.T), flux(self.elem_type), temperature
        )

    def initModel(self, method, material_file):
        super().initModel(method, material_file)

        if method != akantu._static:
            self.model.setTimeStep(0.5 * self.model.getStableTimeStep())


def run_test_generic(self_, method):
    self_.initModel(method, "heat_transfer_input.dat")

    coordinates = self_.mesh.getNodes()
    temperature = self_.model.getTemperature()
    #  set the position of all nodes to the static solution
    self_.setLinearDOF(temperature, coordinates)

    for s in range(0, 100):
        self_.model.solveStep()

    self_.checkAll()
