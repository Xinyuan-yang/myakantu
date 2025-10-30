#!/usr/bin/env python3
__copyright__ = (
    "Copyright (©) 2018-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)"
    "Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)"
)
__license__ = "LGPLv3"


import sys
import akantu
from patch_test_linear_heat_transfer_fixture import TestPatchTestHTMLinear
from patch_test_linear_heat_transfer_fixture import run_test_generic


def run_test(self_):
    run_test_generic(self_, akantu._implicit_dynamic)


def test():
    TestPatchTestHTMLinear.TYPED_TEST(run_test, "Explicit")


if "pytest" not in sys.modules:
    test()
