__copyright__ = (
    "Copyright (©) 2018-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)"
    "Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)"
)
__license__ = "LGPLv3"


import warnings as _aka_warn

import numpy as _aka_np
import scipy.sparse as _aka_sparse

from . import py11_akantu as _py11_akantu

private_keys = set(dir(_py11_akantu)) - set(dir())

for k in private_keys:
    globals()[k] = getattr(_py11_akantu, k)

if _py11_akantu.has_mpi():
    try:
        from mpi4py import MPI  # noqa: F401
    except Exception:
        pass


def initialize(*args, **kwargs):
    raise RuntimeError(
        "No need to call initialize, use parseInput to read an input file"
    )


def finalize(*args, **kwargs):
    _aka_warn.warn("No need to call finalize", DeprecationWarning)


class AkantuSparseMatrix(_aka_sparse.coo_matrix):
    """Wrapper for Akantu sparse matrices that converts to SciPy COO format.

    This class converts Akantu's internal sparse matrix representation to
    SciPy's coordinate (COO) format.

    Args:
        aka_sparse: Akantu sparse matrix AIJ

    Note:
        For symmetric matrices, only the lower triangular part is stored in
        Akantu format. This class automatically expands to full storage by
        mirroring non-diagonal elements.
    """

    def __init__(self, aka_sparse):
        self.aka_sparse = aka_sparse
        matrix_type = self.aka_sparse.getMatrixType()
        sz = self.aka_sparse.size()
        row = self.aka_sparse.getIRN()[:, 0] - 1
        col = self.aka_sparse.getJCN()[:, 0] - 1
        data = self.aka_sparse.getA()[:, 0]

        row = row.copy()
        col = col.copy()
        data = data.copy()

        if matrix_type == _py11_akantu._symmetric:
            non_diags = row != col
            row_sup = col[non_diags]
            col_sup = row[non_diags]
            data_sup = data[non_diags]
            col = _aka_np.concatenate((col, col_sup))
            row = _aka_np.concatenate((row, row_sup))
            data = _aka_np.concatenate((data, data_sup))

        _aka_sparse.coo_matrix.__init__(
            self, (data, (row, col)), shape=(sz, sz), dtype=data.dtype
        )


def registerCustomMaterials(materials):
    """Register custom materials to the MaterialFactory.

    Args: materials (dict): Dictionary mapping material names to material
        classes. Each material class should be callable with (model, id_)
        arguments.

    Example:
        >>> class CustomMaterial(aka.Material):
        ...     def __init__(self, model, id_):
        ...         super().__init__(model, id_)
        >>>
        >>> materials = {"custom_elastic": CustomMaterial}
        >>> registerCustomMaterials(materials)

    """
    mat_factory = _py11_akantu.MaterialFactory.getInstance()

    for name, material in materials.items():
        mat_factory.registerAllocator(
            name, lambda dim, option, model, id_, mat=material: mat(model, id_)
        )


def convertGmshToAkantuMesh(gmsh):  # NOLINT
    """Convert gmsh defined mesh to akantu Mesh.

    Args:
        gmsh: GMSH model object with mesh data

    Returns:
        akantu.Mesh: Converted Akantu mesh object

    Raises:
        ValueError: If mesh is empty or contains unsupported element types
        IndexError: If GMSH element type is not supported
    """
    _msh_to_akantu_element_types = [
        _py11_akantu._not_defined,
        _py11_akantu._segment_2,
        _py11_akantu._triangle_3,
        _py11_akantu._quadrangle_4,
        _py11_akantu._tetrahedron_4,
        _py11_akantu._hexahedron_8,
        _py11_akantu._pentahedron_6,
        _py11_akantu._not_defined,
        _py11_akantu._segment_3,
        _py11_akantu._triangle_6,
        _py11_akantu._not_defined,
        _py11_akantu._tetrahedron_10,
        _py11_akantu._not_defined,
        _py11_akantu._not_defined,
        _py11_akantu._not_defined,
        _py11_akantu._point_1,
        _py11_akantu._quadrangle_8,
        _py11_akantu._hexahedron_20,
        _py11_akantu._pentahedron_15,
    ]

    gmsh_mesh = gmsh.model.mesh
    types = gmsh_mesh.get_element_types()

    if not types:
        raise ValueError("Empty mesh: no element types found")

    dim = -1
    for el_type in types:
        props = gmsh_mesh.get_element_properties(el_type)
        dim = max(dim, props[1])

    if dim < 0:
        raise ValueError("Invalid mesh dimension")

    aka_mesh = _py11_akantu.Mesh(dim)
    accessor = _py11_akantu.MeshAccessor(aka_mesh)

    gmsh_nodes = gmsh_mesh.get_nodes()

    nb_nodes = gmsh_nodes[0].shape[0]
    if nb_nodes == 0:
        raise ValueError("Empty mesh: no nodes found")

    accessor.resizeNodes(nb_nodes)

    aka_nodes = aka_mesh.getNodes()
    gmsh_aka_nodes = {}
    for i in range(nb_nodes):
        aka_nodes[i, :] = gmsh_nodes[1][i * 3 : i * 3 + dim]
        gmsh_aka_nodes[gmsh_nodes[0][i]] = i

    gmsh_aka_elements = {}
    for el_type in types:
        if el_type >= len(_msh_to_akantu_element_types):
            raise IndexError(f"Unsupported GMSH element type: {el_type}")
        aka_type = _msh_to_akantu_element_types[el_type]
        if aka_type == _py11_akantu._not_defined:
            raise ValueError(
                f"GMSH element type {el_type}not supported in Akantu"
            )
        gmsh_aka_elements[aka_type] = {}
        tags, gmsh_conn = gmsh_mesh.get_elements_by_type(el_type)
        aka_mesh.addConnectivityType(aka_type)
        accessor.resizeConnectivity(tags.shape[0], aka_type)
        conn = aka_mesh.getConnectivity(aka_type)

        for element, tag in enumerate(tags):
            gmsh_aka_elements[aka_type][tag] = element
            nnodes_per_el = conn.shape[1]
            for node in range(nnodes_per_el):
                conn[element, node] = gmsh_aka_nodes[
                    gmsh_conn[element * nnodes_per_el + node]
                ]

    for group in gmsh.model.get_physical_groups():
        name = gmsh.model.get_physical_name(*group)
        types = gmsh_mesh.get_element_types(*group)
        aka_group = aka_mesh.createElementGroup(
            name, spatial_dimension=group[0], replace_group=True
        )
        entities = gmsh.model.get_entities_for_physical_group(*group)
        for el_type in types:
            aka_type = _msh_to_akantu_element_types[el_type]
            for entity in entities:
                elements = gmsh_mesh.get_elements_by_type(el_type, tag=entity)
                for el in elements[0]:
                    aka_group.add(
                        _py11_akantu.Element(
                            aka_type,
                            gmsh_aka_elements[aka_type][el],
                            _py11_akantu._not_ghost,
                        ),  # noqa
                        True,
                        False,
                    )
        aka_group.optimize()
    accessor.makeReady()
    return aka_mesh


FromStress = _py11_akantu.FromHigherDim
FromTraction = _py11_akantu.FromSameDim
_py11_akantu.__initialize()

__version__ = _py11_akantu.getVersion()
