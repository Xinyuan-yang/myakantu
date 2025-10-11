#!/usr/bin/env python3
__copyright__ = (
    "Copyright (©) 2019-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)"
    "Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)"
)
__license__ = "LGPLv3"


import akantu as aka


def getStiffnessMatrix(material_file, mesh_file, traction):
    aka.parseInput(material_file)
    spatial_dimension = 2

    # --------------------------------------------------------------------------
    # Initialization
    # --------------------------------------------------------------------------
    mesh = aka.Mesh(spatial_dimension)
    mesh.read(mesh_file)

    model = aka.SolidMechanicsModel(mesh)
    model.initFull(_analysis_method=aka._static)

    model.assembleStiffnessMatrix()
    K = model.getDOFManager().getMatrix("K")
    stiff = aka.AkantuSparseMatrix(K).toarray()

    return stiff


# --------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------
def main():
    mesh_file = "plate.msh"
    material_file = "material.dat"

    traction = 1.0
    mat = getStiffnessMatrix(material_file, mesh_file, traction)
    print(mat)


# --------------------------------------------------------------------------
if __name__ == "__main__":
    main()
