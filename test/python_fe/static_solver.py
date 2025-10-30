"""static_solver.py: static solver for test purposes"""

from . import export
from . import fe

import copy
import scipy.sparse.linalg as spl

# ------------------------------------------------------------------------------
__author__ = "Nicolas Richart"
__credits__ = [
    "Nicolas Richart <nicolas.richart@epfl.ch>",
]
__copyright__ = (
    "Copyright (©) 2016-2021 EPFL (Ecole Polytechnique Fédérale"
    " de Lausanne) Laboratory (LSMS - Laboratoire de Simulation"
    " en Mécanique des Solides)"
)
__license__ = "LGPLv3"
# ------------------------------------------------------------------------------


@export
class StaticSolver(fe.Solver):
    def __init__(self, model):
        self._model = model

    def _assembleResidual(self):
        self._r = self._model.f_ext - self._model.f_int

    def _assembleJacobian(self):
        self._J = copy.copy(self._model.K)
        self._model.applyDirichletBC()

        self._zero_rows(self._J, self._model.blocked)

    def solveStep(self):
        self._assembleJacobian()
        self._assembleResidual()
        x = spl.spsolve(self._J, self._r)
        self._model.u += x
