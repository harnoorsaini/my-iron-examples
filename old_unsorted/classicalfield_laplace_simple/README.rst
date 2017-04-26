=====================================
Laplace - Using OpenCMISS from Python
=====================================

In this tutorial we will walk through how to solve Laplace's equation on a 3D geometry using the Python bindings to OpenCMISS.

-------------
What You Need
-------------

You need to install Python and the OpenCMISS library as described elsewhere.

The Python code given here follows the `Python Laplace example`_ found in the OpenCMISS examples repository.

------------------------------------
Setting Up the Environment Variables
------------------------------------

Environment variables control the programmer setup options. If we are using the bash shell, setting up the environment variables can be done by running the following command from console::

  . ~/.bashrc

And if we are using the c shell::

  . ~/.cshrc

---------------
Getting Started
---------------

In order to use OpenCMISS we have to first import the ``CMISS`` module from the ``opencmiss`` package:

.. literalinclude:: Python/LaplaceExample.py
  :linenos:
  :start-after: # Intialise OpenCMISS start
  :end-before: # Intialise OpenCMISS end

Assuming OpenCMISS has been correctly built with the Python bindings by following the instructions in the `programmer documentation`_, we can now access all the OpenCMISS functions, classes and constants under the ``CMISS`` namespace.

.. _programmer documentation: http://cmiss.bioeng.auckland.ac.nz/OpenCMISS/doc/programmer/

-------------------------------------------
Getting the Computational Nodes Information
-------------------------------------------

OpenCMISS is designed to solve problems on distributed parallel computers. It divides the problem (`decomposes` it) into smaller parts that each run on separate machines/processes, termed `computational nodes`, which solve their part and intercommunicate to solve the whole problem. We can get the computational nodes number as below:

.. literalinclude:: Python/LaplaceExample.py
  :linenos:
  :start-after: #DOC-START computationalNode
  :end-before: #DOC-END computationalNode

----------------------------
Creating a Coordinate System
----------------------------

First we construct a `coordinate system` that will be used to describe the geometry in our problem. The 3D geometry will exist in a 3D space, so we need a 3D coordinate system.

When creating an object in OpenCMISS there are at least three steps. First we initialise the object:

.. literalinclude:: Python/LaplaceExample.py
  :linenos:
  :start-after: #DOC-START coordinate1
  :end-before: #DOC-END coordinate1

This creates a thin wrapper that just points to the actual coordinate system used internally by OpenCMISS, and initially the pointer is null. Trying to use the coordinate system now would raise an exception. To actually construct a coordinate system we call the ``CreateStart`` method and pass it a user number. The user number must be unique between objects of the same type and can be used to identify the coordinate system later. Most OpenCMISS classes require a user number when creating them, and many also require additional parameters to the ``CreateStart`` method:

.. literalinclude:: Python/LaplaceExample.py
  :linenos:
  :start-after: #DOC-START coordinate2
  :end-before: #DOC-END coordinate2

We can now optionally set any properties on the object. We will set the dimension so that the coordinate system is 3D:

.. literalinclude:: Python/LaplaceExample.py
  :linenos:
  :start-after: #DOC-START coordinate3
  :end-before: #DOC-END coordinate3

And finally, we finish creating the coordinate system:

.. literalinclude:: Python/LaplaceExample.py
  :linenos:
  :start-after: #DOC-START coordinate4
  :end-before: #DOC-END coordinate4

-----------------
Creating a Region
-----------------

Next we create a `region` that our fields will be defined on and tell it to use the 3D coordinate system we created previously. The ``CreateStart`` method for a region requires another region as a parameter. We use the world region that is created by default so that our region is a sub-region of the world region. We can also give the region a label:

.. literalinclude:: Python/LaplaceExample.py
  :linenos:
  :start-after: #DOC-START region
  :end-before: #DOC-END region

----------------
Creating a Basis
----------------

The finite element description of our fields requires a `basis function` to interpolate field values over elements, so we create a 3D basis with linear Lagrange interpolation in all three Xi directions:

.. literalinclude:: Python/LaplaceExample.py
  :linenos:
  :start-after: #DOC-START basis
  :end-before: #DOC-END basis

When we set the basis type we select a value from the ``BasisTypes`` enum.

--------------------------
Creating a Decomposed Mesh
--------------------------

In order to define a simple 3D geometry for our problem we can use one of OpenCMISS's inbuilt generated meshes. We will create a 3D, cuboid mesh with 5 elements in the X, Y and Z directions and tell it to use the basis we created previously:

.. literalinclude:: Python/LaplaceExample.py
  :linenos:
  :start-after: #DOC-START generated mesh
  :end-before: #DOC-END generated mesh

When setting the ``basis`` property, we assign a list of bases as we might want to construct a mesh with multiple components using different interpolation schemes. In this example we only have one mesh component.

The generated mesh is not itself a mesh, but is used to create a mesh. We construct the `mesh` when we call the ``CreateFinish`` method of the generated mesh and pass in the mesh to generate:

.. literalinclude:: Python/LaplaceExample.py
  :linenos:
  :start-after: #DOC-START mesh
  :end-before: #DOC-END mesh

Here we have initialised a mesh but not called ``CreateStart`` or ``CreateFinish``, instead the mesh creation is done when finishing the creation of the generated mesh.

The next step in this example is to decompose the mesh for the number of computation nodes in use i.e. creating a decomposition with that number of sub-domains:

.. literalinclude:: Python/LaplaceExample.py
  :linenos:
  :start-after: #DOC-START decomposition
  :end-before: #DOC-END decomposition

Note that even when we have just one computational node, OpenCMISS still needs to work with a decomposed mesh, which will have one domain.

-----------------
Defining Geometry
-----------------

Now that we have a decomposed mesh, we can begin defining the `fields` we need on it. First we will create a geometric field to define our problem geometry:

.. literalinclude:: Python/LaplaceExample.py
  :linenos:
  :start-after: #DOC-START geometry
  :end-before: #DOC-END geometry

The call to the ``ComponentMeshComponentSet`` method is not actually required here as all field components will default to use the first mesh component that we created with ``generatedMesh.CreateStart`` method, but if we have defined a mesh that has multiple components (that use different interpolation schemes) then different field components can use different mesh components. For example, in a finite elasticity problem we could define our geometry using quadratic Lagrange interpolation, and the hydrostatic pressure using linear Lagrange interpolation.

We have created a field but all the field component values are currently set to zero. We can define the geometry using the generated mesh we created earlier:

.. literalinclude:: Python/LaplaceExample.py
  :linenos:
  :start-after: #DOC-START geometry1
  :end-before: #DOC-END geometry2

--------------------
Setting up Equations
--------------------

Now we have a geometric field we can construct an `equations set`. This defines the set of equations that we wish to solve in our problem on this region. The specific equation set we are solving is defined by ``equationsSetSpecification`` list which is the fourth parameter to the ``CreateStart`` method. The first, second and third parameters in the list are the equations set class, type and subtype respectively. In this example we are solving the standard Laplace equation which is a member of the classical field equations set class and the Laplace equation type. When we create an equations set we also have to create an equations set field, however, this is only used to identify multiple equations sets of the same type on a region so we will not use it:

.. literalinclude:: Python/LaplaceExample.py
  :linenos:
  :start-after: #DOC-START equationset
  :end-before: #DOC-END equationset

Now we use our equations set to create a dependent field. This stores the solution to our equations:

.. literalinclude:: Python/LaplaceExample.py
  :linenos:
  :start-after: #DOC-START dependent
  :end-before: #DOC-END dependent

We haven't used the ``Field.CreateStart`` method to construct the dependent field but have had it automatically constructed by the equations set.

We can initialise our solution with a value we think will be close to the final solution. A field in OpenCMISS can contain multiple `field variables`, and each field variable can have multiple `components`. For the standard Laplace equation, the dependent field has ``U`` (standard variable type i.e., u) and ``DELUDELN`` (normal derivative variable type i.e., du/dn) variables which both have one component. Field variables can also have different field `parameter sets`, for example we can store values at a previous time step in dynamic problems. In this example we are only interested in the ``VALUES`` parameter set:

.. literalinclude:: Python/LaplaceExample.py
  :linenos:
  :start-after: #DOC-START dependent1
  :end-before: #DOC-END dependent1

Once the equations set is defined, we create the `equations` that use our fields to construct equations matrices and vectors. We will use sparse matrices to store the equations and disable output when assembling the equations:

.. literalinclude:: Python/LaplaceExample.py
  :linenos:
  :start-after: #DOC-START set
  :end-before: #DOC-END set

--------------------
Defining the Problem
--------------------

Now that we have defined all the equations we will need we can create our `problem` to solve. We create a standard Laplace problem, which is a member of the classical field problem class and Laplace equation problem type:

.. literalinclude:: Python/LaplaceExample.py
  :linenos:
  :start-after: #DOC-START problem
  :end-before: #DOC-END problem

The problem type defines a `control loop` structure that is used when solving the problem. We may have multiple control loops with nested sub-loops, and control loops can have different types, for example load incremented loops or time loops for dynamic problems. In this example a simple, single iteration loop is created without any sub-loops. If we wanted to access the control loop and modify it we would use the ``problem.ControlLoopGet`` method before finishing the creation of the control loops, but we will just leave it with the default configuration:

.. literalinclude:: Python/LaplaceExample.py
  :linenos:
  :start-after: #DOC-START loops
  :end-before: #DOC-END loops

-------------------
Configuring Solvers
-------------------

After defining the problem structure we can create the `solvers` that will be run to actually solve our problem. The problem type defines the solvers to be set up so we call ``problem.SolversCreateStart`` to create the solvers and then we can access the solvers to modify their properties:

.. literalinclude:: Python/LaplaceExample.py
  :linenos:
  :start-after: #DOC-START problem solver
  :end-before: #DOC-END problem solver

Note that we initialised a solver, created it with the call to ``SolversCreateStart`` and then we obtained it with the call to ``SolverGet``. If we look at the help for the ``SolverGet`` method we see it takes three parameters:

controlLoopIdentifiers:
    A list of integers used to identify the control loop to get a solver for.
    This always starts with the root control loop, given by ``CMISS.ControlLoopIdentifiers.NODE``.
    In this example we only have the one control loop and no sub-loops.

solverIndex:
    The index of the solver to get, as a control loop may have multiple solvers.
    In this case there is only one solver in our root control loop.

solver:
    An initialised solver object that hasn't been created yet, and on return
    it will be the solver that we asked for.

Once we've obtained the solver we then set various properties before finishing the creation of all the problem solvers.

After defining our solver we can create the equations for the solver to solve by adding our equations sets to the solver equations. In this example we have just one equations set to add but for coupled problems we may have multiple equations sets in the solver equations. We also tell OpenCMISS to use sparse matrices to store our solver equations:

.. literalinclude:: Python/LaplaceExample.py
  :linenos:
  :start-after: #DOC-START solver equations
  :end-before: #DOC-END solver equations

---------------------------
Setting Boundary Conditions
---------------------------

The final step in configuring the problem is to define the boundary conditions to be satisfied. We will set the dependent field value at the first node to be zero, and at the last node to be 1.0. These nodes will correspond to opposite corners in our geometry. Because OpenCMISS can solve our problem on multiple computational nodes where each computational node does not necessarily know about all nodes in our mesh, we must first check that the node we are setting the boundary condition at is in our computational node domain:

.. literalinclude:: Python/LaplaceExample.py
  :linenos:
  :start-after: #DOC-START BC
  :end-before: #DOC-END BC

In parallel, the exact same code is running on each computational node. When setting a boundary condition at a node we can use either the ``AddNode`` method or the ``SetNode`` method. Using ``AddNode`` will add the value we provide to the current field value and set this as the boundary condition value, but here we want to directly specify the value so we use the ``SetNode`` method.

The arguments to the ``SetNode`` method are the field, field variable type, node version number, node derivative number, node user number, field component number, boundary condition type and boundary condition value. The version and derivative numbers are one as we aren't using versions and we are setting field values rather than derivative values. We can also only set derivative boundary conditions when using a Hermite basis type. There are a wide number of boundary condition types that can be set but many are only available for certain equation set types and in this example we simply want to fix the field value.

When ``solverEquations.BoundaryConditionsCreateFinish()`` is called OpenCMISS will construct the solver matrices and vectors.

Since the Laplace equation is the steady-state heat equation, one physical interpretation of this problem is as follows: fix the temperature on the boundary of the domain according to the given specification of the boundary condition. Allow heat to flow until a stationary state is reached in which the temperature at each point on the domain doesn't change anymore. The temperature distribution in the interior will then be given by the solution to the corresponding Dirichlet problem.

-------
Solving
-------

After our problem solver equations have been fully defined we are now ready to solve our problem. When we call the ``Solve`` method of the problem it will loop over the control loops and control loop solvers to solve our problem:

.. literalinclude:: Python/LaplaceExample.py
  :linenos:
  :start-after: #DOC-START Solve
  :end-before: #DOC-END Solve

----------------------
Exporting the Solution
----------------------

Once the problem has been solved, the dependent field contains the solution to our problem. We can then export the dependent and geometric fields to a FieldML file so that we can visualise the solution using cmgui. We will export the geometric and dependent field values to a ``LaplaceExample.xml`` file. Separate plain text data files will also be created:

.. literalinclude:: Python/LaplaceExample.py
  :linenos:
  :start-after: #DOC-START Export
  :end-before: #DOC-END Export

===========================
Running the Script in Linux
===========================

The Python code given here follows the `Python Laplace example`_ found in the OpenCMISS examples repository. Now we would like to run this script to check that we get the correct output from OpenCMISS. In each case you must change directory to where you downloaded the Python script and data file, then run the script with ``python``. The exact names of the directories may not match what is on your own computer, so you will need to change them as appropriate. Follow the instructions that are applicable for your platform.

.. _Python Laplace example: https://github.com/OpenCMISS-Examples/laplace/blob/master/Python/LaplaceExample.py
                            
------
Serial
------

From console::

  python LaplaceExample.py

--------
Parallel
--------

From console::

  cd ~/.
  touch .mpd.conf

Edit .mpd.conf and add::

  MPD_SECRETWORD=<your secret word>

DO NOT use your password for your secret word. Now::

  chmod 600 .mpd.conf

Now edit a file called hostfile.list and add the following lines (one line per processor)::

  hpc3.bioeng.auckland.ac.nz

Finally add the following to your login script. For bash::

  export MP_HOSTFILE=<path to your hostfile.list file> 

For cshell::

  setenv MP_HOSTFILE <path to your hostfile.list file> 

Run a MPI example using mpiexec::

  cd ${OPENCMISS_ROOT}/ClassicalField/Laplace/Laplace/Python
  mpd &
  mpiexec.mpd –n 2 python LaplaceExample.py

To see what is happening get another terminal up on hpc3 and type top.

======
Output
======
  
If OpenCMISS is installed and running correctly then you should see the solver matrix as output in the console window similar to::

      [0.454768E+00  0.453061E+00  0.474333E+00  0.479934E+00  0.482290E+00  0.316870E+00  0.409378E+00  0.462585E+00
       0.472722E+00  0.480318E+00  0.482162E+00  0.402008E+00  0.428823E+00  0.460520E+00  0.473907E+00  0.480185E+00
       0.482271E+00  0.433056E+00  0.443886E+00  0.462155E+00  0.474202E+00  0.480324E+00  0.482290E+00  0.445737E+00
       0.451788E+00  0.464094E+00  0.474444E+00  0.480395E+00  0.482433E+00  0.449290E+00  0.454204E+00  0.464855E+00
       0.474531E+00  0.480518E+00  0.482091E+00  0.515997E+00  0.439936E+00  0.470055E+00  0.478714E+00  0.484677E+00
       0.486369E+00  0.450287E+00  0.454469E+00  0.467957E+00  0.479295E+00  0.484600E+00  0.486426E+00  0.450036E+00
       0.457122E+00  0.469658E+00  0.479334E+00  0.484718E+00  0.486441E+00  0.455945E+00  0.460982E+00  0.471032E+00
       0.479634E+00  0.484782E+00  0.486447E+00  0.460150E+00  0.463890E+00  0.472092E+00  0.479883E+00  0.484851E+00
       0.486268E+00  0.461582E+00  0.464923E+00  0.472480E+00  0.480031E+00  0.484658E+00  0.487089E+00  0.471228E+00
       0.488941E+00  0.487827E+00  0.493293E+00  0.496779E+00  0.498081E+00  0.486044E+00  0.485514E+00  0.488588E+00
       0.493186E+00  0.496810E+00  0.498062E+00  0.484843E+00  0.485417E+00  0.488683E+00  0.493308E+00  0.496818E+00
       0.498090E+00  0.483959E+00  0.485250E+00  0.488910E+00  0.493397E+00  0.496893E+00  0.498314E+00  0.483784E+00
       0.485276E+00  0.489069E+00  0.493481E+00  0.496932E+00  0.498933E+00  0.483796E+00  0.485308E+00  0.489169E+00
       0.493323E+00  0.497789E+00  0.495580E+00  0.504420E+00  0.502211E+00  0.506677E+00  0.510831E+00  0.514692E+00
       0.516204E+00  0.501067E+00  0.503068E+00  0.506519E+00  0.510931E+00  0.514724E+00  0.516216E+00  0.501686E+00
       0.503107E+00  0.506603E+00  0.511090E+00  0.514750E+00  0.516041E+00  0.501910E+00  0.503182E+00  0.506692E+00
       0.511317E+00  0.514583E+00  0.515157E+00  0.501938E+00  0.503190E+00  0.506814E+00  0.511412E+00  0.514486E+00
       0.513956E+00  0.501919E+00  0.503221E+00  0.506707E+00  0.512173E+00  0.511059E+00  0.528772E+00  0.512911E+00
       0.515342E+00  0.519969E+00  0.527520E+00  0.535077E+00  0.538418E+00  0.513732E+00  0.515149E+00  0.520117E+00
       0.527908E+00  0.536110E+00  0.539850E+00  0.513553E+00  0.515218E+00  0.520366E+00  0.528968E+00  0.539018E+00
       0.544055E+00  0.513559E+00  0.515282E+00  0.520666E+00  0.530342E+00  0.542878E+00  0.549964E+00  0.513574E+00
       0.515400E+00  0.520705E+00  0.532043E+00  0.545531E+00  0.549713E+00  0.513631E+00  0.515323E+00  0.521286E+00
       0.529945E+00  0.560064E+00  0.484003E+00  0.517909E+00  0.519482E+00  0.525469E+00  0.535145E+00  0.545796E+00
       0.550710E+00  0.517567E+00  0.519605E+00  0.525556E+00  0.535906E+00  0.548212E+00  0.554263E+00  0.517710E+00
       0.519676E+00  0.525798E+00  0.537845E+00  0.556114E+00  0.566944E+00  0.517729E+00  0.519815E+00  0.526093E+00
       0.539480E+00  0.571177E+00  0.597992E+00  0.517838E+00  0.519682E+00  0.527278E+00  0.537415E+00  0.590622E+00
       0.683130E+00  0.517710E+00  0.520066E+00  0.525667E+00  0.546939E+00  0.545232E+00]

If this is your first use of OpenCMISS, *congratulations on getting this far!*
