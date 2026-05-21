.. _sect-io:

Input/Output
============

Input file
----------

The text input file of a simulation should be precised using the method
:cpp:func:`initialize <akantu::initialize>` which will instantiate the static
:cpp:class:`Parser <akantu::Parser>` object of ``Akantu``. This section
explains how to manipulate at :cpp:class:`Parser <akantu::Parser>` objects to
input data in ``Akantu``.

Akantu Parser
~~~~~~~~~~~~~

``Akantu`` file parser has a tree organization.

- :cpp:class:`Parser <akantu::Parser>`, the root of the tree, can be accessed
  using::

     auto & parser = getStaticParser();

- :cpp:class:`ParserSection <akantu::ParserSection>`, branch of the tree,
  contains map a of sub-sections (:cpp:enum:`ParserType
  <akantu::ParserType>`, :cpp:class:`ParserSection <akantu::ParserSection>`)
  and a :cpp:class:`ParserSection * <akantu::ParserSection>` pointing to the
  parent section. The user section of the input file can directly be accessed
  by::

     const auto & usersect = getUserParser();

- :cpp:class:`ParserParameter <akantu::ParserParameter>`, the leaf of the tree,
  carries data of the input file which can be cast to the correct type with the
  assignment operator::

     Real mass = usersect.getParameter("mass");

  or used directly within an expression

Grammar
~~~~~~~

The structure of text input files consists of different sections
containing a list of parameters. As example, the file parsed in the
previous section will look like::

  user parameters [
    mass = 10.5
  ]

Basically every standard arithmetic operations can be used inside of input files
as well as the constant ``pi`` and ``e`` and the exponent operator ``^``.
Operations between :cpp:class:`ParserParameter <akantu::ParserParameter>` are
also possible with the convention that only parameters of the current and the
parent sections are available. :cpp:class:`Vector <akantu::Vector>` and
:cpp:class:`Matrix <akantu::Matrix>` can also be read according to the ``NumPy``
:cite:`numpy` writing convention (a.e. cauchy_stress_tensor = [[:math:`\sigma_{xx}`,
:math:`\sigma_{xy}`],[:math:`\sigma_{yx}`,\ :math:`\sigma_{yy}`]]). An example
illustrating how to parse the following input file can be found in
``example\io\parser\example_parser.cc``::

  user parameters [
      spatial_dimension = 2
      mesh_file     = swiss_cheese.msh
      inner_holes   = holes
      outter_crust  = crust
      lactostatic_p = 30e3
      stress        = [[lactostatic_p, 0            ],
                       [0,             lactostatic_p]]
      max_nb_iterations = 100
      precision     = 1e-9
  ]

.. _sect-io-solver:

Model solver configuration
~~~~~~~~~~~~~~~~~~~~~~~~~~

Solver settings are read from a ``model`` section whose name must match the
model type (e.g. ``solid_mechanics_model``).  An optional ``name`` parameter
can be used to disambiguate when several models of the same type coexist::

  model solid_mechanics_model [
    dof_manager_type = petsc   # "default" (MUMPS/Eigen) or "petsc"
    default_solver   = static  # name of the time_step_solver to activate

    time_step_solver static [
      ...
    ]
  ]

.. rubric:: Model-level parameters

.. list-table::
   :header-rows: 1
   :widths: 25 15 60

   * - Parameter
     - Default
     - Description
   * - ``dof_manager_type``
     - ``default``
     - DOF manager backend: ``default`` (MUMPS/Eigen) or ``petsc``
   * - ``default_solver``
     - *(first defined)*
     - Name of the ``time_step_solver`` block to use as the active solver

.. rubric:: ``time_step_solver <type>`` sub-section

``<type>`` is one of ``static``, ``dynamic``, or ``dynamic_lumped``.

.. list-table::
   :header-rows: 1
   :widths: 25 15 60

   * - Parameter
     - Default
     - Description
   * - ``name``
     - *type string*
     - Optional alias used to reference this solver (e.g. in ``default_solver``)
   * - ``sparse_solver_type``
     - ``auto``
     - Sparse-matrix backend: ``mumps``, ``eigen``, ``petsc``, or ``auto``

A ``time_step_solver`` block may contain at most one ``non_linear_solver``
sub-section and any number of ``integration_scheme`` sub-sections.

.. rubric:: ``non_linear_solver <type>`` sub-section

``<type>`` is one of ``linear``, ``newton_raphson``,
``newton_raphson_modified``, ``lumped``, ``petsc_snes``, or ``petsc_tao``.

**Newton-Raphson types** (``newton_raphson``, ``newton_raphson_modified``):

.. list-table::
   :header-rows: 1
   :widths: 30 15 55

   * - Parameter
     - Default
     - Description
   * - ``max_iterations``
     - 10
     - Maximum number of Newton iterations
   * - ``threshold``
     - 1e-10
     - Convergence threshold
   * - ``convergence_type``
     - ``solution``
     - Convergence criterion: ``residual``, ``solution``, or
       ``residual_mass_wgh``
   * - ``force_linear_recompute``
     - ``true``
     - Reassemble the Jacobian at every iteration
   * - ``verbose_iteration``
     - ``false``
     - Print convergence information at each iteration

**PETSc types** (``petsc_snes``, ``petsc_tao``):

All parameters are forwarded verbatim to PETSc's option database
(prefixed with ``-``), so any PETSc option accepted by ``SNESSetFromOptions``
or ``TaoSetFromOptions`` is valid.  Command-line PETSc options take
precedence and are applied on top of the file settings.

The Akantu-side names ``max_iterations``, ``threshold``, and
``convergence_type`` are **not** accepted for PETSc solvers; use the
corresponding PETSc options instead:

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - PETSc SNES option
     - Meaning
   * - ``snes_max_it``
     - Maximum number of SNES iterations
   * - ``snes_atol``
     - Absolute residual tolerance
   * - ``snes_rtol``
     - Relative residual tolerance
   * - ``snes_stol``
     - Step-length tolerance
   * - ``snes_type``
     - SNES algorithm — see the
       `SNES manual <https://petsc.org/release/manual/snes/>`_ and
       `SNESType <https://petsc.org/release/manualpages/SNES/SNESType/>`_
       (e.g. ``newtonls``, ``ncg``, ``nrichardson``)
   * - ``ksp_type``
     - Linear solver inside SNES — see the
       `KSP manual <https://petsc.org/release/manual/ksp/>`_ and
       `KSPType <https://petsc.org/release/manualpages/KSP/KSPType/>`_
       (e.g. ``cg``, ``gmres``, ``preonly``)
   * - ``pc_type``
     - Preconditioner — see the
       `KSP manual <https://petsc.org/release/manual/ksp/>`_ and
       `PCType <https://petsc.org/release/manualpages/PC/PCType/>`_
       (e.g. ``hypre``, ``ilu``, ``lu``)

For TAO solvers (``petsc_tao``) the equivalent options are prefixed with
``tao_`` — see the
`TAO manual <https://petsc.org/release/manual/tao/>`_ and
`TaoType <https://petsc.org/release/manualpages/Tao/TaoType/>`_
for the full list of algorithms.

Example for a PETSc Newton-Krylov solve::

  model solid_mechanics_model [
    dof_manager_type = petsc
    default_solver   = static

    time_step_solver static [
      non_linear_solver petsc_snes [
        snes_max_it = 500
        snes_atol   = 1e-8
        snes_rtol   = 1e-8
        pc_type     = hypre
        ksp_type    = cg
      ]
    ]
  ]

.. rubric:: ``integration_scheme <type>`` sub-section

Used with dynamic solvers to specify the time-integration scheme for each DOF
type.  ``<type>`` is one of ``pseudo_time``, ``forward_euler``,
``backward_euler``, ``central_difference``, ``trapezoidal_rule_1``,
``trapezoidal_rule_2``, ``newmark_beta``, ``generalized_trapezoidal``, etc.

.. list-table::
   :header-rows: 1
   :widths: 25 15 60

   * - Parameter
     - Default
     - Description
   * - ``name``
     - *(required)*
     - DOF identifier this scheme applies to (e.g. ``displacement``)
   * - ``solution_type``
     - *(model default)*
     - Correction variable: ``displacement``, ``velocity``, or ``acceleration``

.. _sect-io-material:

Material section
~~~~~~~~~~~~~~~~

Material characteristics (constitutive behaviour and properties) are specified
with ``material`` sub-sections.  They can be placed either inside the
``model`` block (preferred) or at the top level of the input file::

  model solid_mechanics_model [
    material constitutive_law [
      name = value
      rho  = value
      ...
    ]
  ]

When no ``model`` block is present the parser falls back to reading
``material`` sections from the global scope::

  material constitutive_law [
    name = value
    rho  = value
    ...
  ]

In both cases *constitutive_law* is the name of the adopted constitutive law,
followed by its properties listed one per line.  Some constitutive laws accept
an *optional flavor*.  More information can be found in
:ref:`sect-smm-cl`.

Output data
-----------

Generic data
~~~~~~~~~~~~

In this section, we address ways to get the internal data in human-readable
formats. The models in ``Akantu`` handle data associated to the mesh, but this
data can be split into several :cpp:class:`Arrays <akantu::Array>`. For example,
the data stored per element type in a :cpp:class:`ElementTypeMapArray
<akantu::ElementTypeMapArray>` is composed of as many :cpp:class:`Arrays
<akantu::Array>` as types in the mesh.

In order to get this data in a visualization software, the models contain a
object to dump ``VTK`` files. These files can be visualized in software such
as ``ParaView`` :cite:`paraview`, ``ViSit`` :cite:`visit` or ``Mayavi``
:cite:`mayavi`.

The internal dumper of the model can be configured to specify which data
fields are to be written. This is done with the  :cpp:func:`addDumpField <akantu::Model::addDumpField>` method. By default all the
files are generated in a folder called ``paraview/``::

  model.setBaseName("output"); // prefix for all generated files
  model.addDumpField("displacement"); model.addDumpField("stress"); ...
  model.dump()

The fields are dumped with the number of components of the memory. For example,
in 2D, the memory has :cpp:class:`Vectors <akantu::Vector>` of 2 components, or
the :math:`2^{nd}` order tensors with :math:`2\times2` components. This memory
can be dealt with :cpp:func:`addDumpFieldVector
<akantu::Model::addDumpFieldVector>` which always dumps :cpp:class:`Vectors
<akantu::Vector>` with 3 components or :cpp:func:`addDumpFieldTensor
<akantu::Model::addDumpFieldTensor>` which dumps :math:`2^{nd}` order tensors
with :math:`3\times3` components respectively. The routines :cpp:func:`addDumpFieldVector <akantu::Model::addDumpFieldVector>` and
:cpp:func:`addDumpFieldTensor <akantu::Model::addDumpFieldTensor>` were
introduced because of ``ParaView`` which mostly manipulate 3D data.

Those fields which are stored by quadrature point are modified to be seen in the
``VTK`` file as elemental data. To do this, the default is to average the
values of all the quadrature points.

The list of fields depends on the models (for :cpp:class:`SolidMechanicsModel
<akantu::SolidMechanicsModel>` see table :ref:`tab-io-smm-field-list`.

.. _tab-io-smm-field-list:

.. table:: List of dumpable fields for :cpp:class:`SolidMechanicsModel <akantu::SolidMechanicsModel>`.

      ====================== ================ =================
      key                    type             support
      ====================== ================ =================
      displacement           ``Vector<Real>`` nodes
      mass                   ``Vector<Real>`` nodes
      velocity               ``Vector<Real>`` nodes
      acceleration           ``Vector<Real>`` nodes
      force                  ``Vector<Real>`` nodes
      residual               ``Vector<Real>`` nodes
      increment              ``Vector<Real>`` nodes
      blocked_dofs           ``Vector<bool>`` nodes
      partitions             ``Real``         elements
      material_index         variable         elements
      strain                 ``Matrix<Real>`` quadrature points
      Green strain           ``Matrix<Real>`` quadrature points
      principal strain       ``Vector<Real>`` quadrature points
      principal Green strain ``Vector<Real>`` quadrature points
      grad_u                 ``Matrix<Real>`` quadrature points
      stress                 ``Matrix<Real>`` quadrature points
      Von Mises stress       ``Real``         quadrature points
      material_index         variable         quadrature points
      ====================== ================ =================

Cohesive elements’ data
~~~~~~~~~~~~~~~~~~~~~~~

Cohesive elements and their relative data can be easily dumped thanks to
a specific dumper contained in :cpp:class:`SolidMechanicsModelCohesive <akantu::SolidMechanicsModelCohesive>`. In order to use it, one has
just to add the string ``cohesive elements`` when calling each method already
illustrated. Here is an example on how to dump displacement and damage::

  model.addDumpFieldVectorToDumper("cohesive elements", "displacement");
  model.addDumpFieldToDumper("cohesive elements", "damage");
  model.dump("cohesive elements");


Fragmentation data
^^^^^^^^^^^^^^^^^^

Whenever the :cpp:class:`SolidMechanicsModelCohesive
<akantu::SolidMechanicsModelCohesive>` is used, it is possible to dump
additional data about the fragments that get formed in the simulation both in
serial and parallel. This task is carried out by the
:cpp:class:`FragmentManager <akantu::FragmentManager>` class, that takes care of
computing the following quantities for each fragment:

-  index;

-  mass;

-  moments of inertia;

-  velocity;

-  number of elements.

These computations can be realized at once by calling the function
:cpp:class:`computeAllData <akantu::FragmentManager::computeAllData>`, or
individually by calling the other public functions of the class. The data can be
dumped to be visualized in ``ParaView``, or can be accessed within the
simulation. An example of usage is:

At the end of this example the velocities of the fragments are accessed with a
reference to a :cpp:class:`const Array\<Real\> <akantu::Array>`. The size of this
array is the number of fragments, and its number of components is the spatial
dimension in this case.

Advanced dumping
~~~~~~~~~~~~~~~~

Arbitrary fields
^^^^^^^^^^^^^^^^

In addition to the predetermined fields from the models and materials,
the user can add any data to a dumper as long as the support is the
same. That is to say data that have the size of the full mesh on if the
dumper is dumping the mesh, or of the size of an element group if it is
a filtered dumper.

For this the easiest is to use the “external” fields register functions

The simple case force nodal and elemental data are to pass directly the
data container itself if it as the good size.

-  For nodal fields:

   It is assumed that the array as the same size as the number of nodes
   in the mesh

-  For elemental fields:

   It is assumed that the arrays in the map have the same sizes as the element
   numbers in the mesh for element types of dimension ``spatial_dimension``.

If some changes have to be applied on the data as for example a padding for
``ParaView`` vectors, this can be done by using the field interface.

All these functions use the default dumper registered in the mesh but also have
the ``ToDumper`` variation with the dumper name specified. For example:

An example of code presenting this interface is present in the
``examples/io/dumper``. This interface is part of the  :cpp:class:`Dumpable
<akantu::Dumpable>` class from which the :cpp:class:`Mesh <akantu::Mesh>`
inherits.

Creating a new dumper
^^^^^^^^^^^^^^^^^^^^^

You can also create you own dumpers, ``Akantu`` uses a third-party library in
order to write the output files, ``IOHelper``. ``Akantu`` supports the
``ParaView`` format and a Text format defined by ``IOHelper``.

This two files format are handled by the classes :cpp:class:`DumperParaview
<akantu::DumperParaview>` and :cpp:class:`DumperText <akantu::DumperText>`.

In order to use them you can instantiate on of this object in your code. This
dumper have a simple interface. You can register a mesh :cpp:func:`registerMesh
<akantu::DumperIOHelper::registerMesh>`, :cpp:func:`registerFilteredMesh
<akantu::DumperIOHelper::registerFilteredMesh>` or a field,
:cpp:class:`registerField <akantu::DumperIOHelper::registerField>`.

An example of code presenting this low level interface is present in the
``examples/io/dumper``. The different types of :cpp:class:`Field
<akantu::Field>` that can be created are present in the source folder
``src/io/dumper``.
