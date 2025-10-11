#!/usr/bin/env python3
__copyright__ = (
    "Copyright (©) 2020-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)"
    "Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)"
)
__license__ = "LGPLv3"


import math


class SinBoundary:
    def __init__(self, amplitude, phase):
        self.amplitude = amplitude
        self.phase = phase

    def compute(self, disp, coord, flags):
        disp[1] = -self.amplitude * math.sin(self.phase * coord[1])
        flags[1] = True
