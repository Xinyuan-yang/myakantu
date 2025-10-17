__copyright__ = (
    "Copyright (©) 2016-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)"
    "Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)"
)
__license__ = "LGPLv3"


import gdb

# Load the pretty-printers.
from .printers import register_akantu_printers

register_akantu_printers(gdb.current_objfile())
