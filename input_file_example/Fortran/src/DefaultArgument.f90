
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! SETTING DEFAULT VALUES OF DATA STRUCTURES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!! A lot of parameters have been left blank " " for  future use !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! default parameters for the solver block
  all_Solver(1:NUmberOfSolver)%NewtonRelativeTolerance(1)                          = "0"
  all_Solver(1:NUmberOfSolver)%NewtonAbsoluteTolerance(1)                          = "0"
  all_Solver(1:NUmberOfSolver)%IterativeSolverAbsoluteTolerance(1)                 = "0"
  all_Solver(1:NUmberOfSolver)%IterativeSolverRelativeTolerance(1)                 = "0"
  all_Solver(1:NUmberOfSolver)%IterativeSolverAbsoluteTolerance(1)                 = "0"
  all_Solver(1:NUmberOfSolver)%JacobianType(1)                                     = "SOLVER_NEWTON_JACOBIAN_EQUATIONS_CALCULATED"
  all_Solver(1:NUmberOfSolver)%LibraryType(1)                                      = "SOLVER_MUMPS_LIBRARY"
  all_Solver(1:NUmberOfSolver)%JacobianType(1)                                     = "SOLVER_ITERATIVE_NO_PRECONDITIONER"
  all_Solver(1:NUmberOfSolver)%LinearSolverType(1)                                 = " "
  all_Solver(1:NUmberOfSolver)%NonLinearSolverType(1)                                 = " "
  all_Solver(1:NUmberOfSolver)%DirectSolverType(1)                                 = " "
  all_Solver(1:NUmberOfSolver)%IterativeSolverType(1)                              = " "
  all_Solver(1:NUmberOfSolver)%IterativeSolverMaximumIterations(1)                 = " "
  all_Solver(1:NUmberOfSolver)%OutputTypeSet(1)                                    = " "

   ! default parameters for the  basis block

  all_Basis(1: NumberOfBasis)%BasisId(1)                                           = " "
  all_Basis(1: NumberOfBasis)%BasisInterpolationType(1)                            = " "
  all_Basis(1: NumberOfBasis)%BasisNumberOfGaussPoints(1)                          = " "
  all_Basis(1: NumberOfBasis)%BasisNumberOfGaussPoints(2)                          = " "
  all_Basis(1: NumberOfBasis)%BasisNumberOfGaussPoints(3)                          = " "

   ! default parameters for the problem block

  all_Problem(1: NUmberOfProblem)%ProblemId(1)                                     = " "
  all_Problem(1: NUmberOfProblem)%ProblemClass(1)                                  = " "
  all_Problem(1: NUmberOfProblem)%ProblemType(1)                                   = " "
  all_Problem(1: NUmberOfProblem)%ProblemSubType(1)                                = " "

   ! default type for the coordinate system block

  all_CoordinateSystem(1:NumberofCoordinateSystem)%CoordinateSystemId(1)           = ""
  all_CoordinateSystem(1:NumberofCoordinateSystem)%CoordinateSystemType(1)         = "COORDINATE_RECTANGULAR_CARTESIAN_TYPE"
  all_CoordinateSystem(1:NumberofCoordinateSystem)%CoordinateSystemDimension(1)    = ""

   ! default type for the coordinate system block

  all_ControlLoop(1:NumberofControlLoop)%ControlLoopId(1)                          = " "
  all_ControlLoop(1:NumberofControlLoop)%ControlLoopType(1)                        = " "
  all_ControlLoop(1:NumberofControlLoop)%ControlLoopTimeIncrement(1)               = " "
  all_ControlLoop(1:NumberofControlLoop)%ControlLoopTimeIncrement(2)               = " "
  all_ControlLoop(1:NumberofControlLoop)%ControlLoopTimeIncrement(3)               = " "
  all_ControlLoop(1:NumberofControlLoop)%ControlLoopLoadIncrement(1)               = " "

  !! default type for the GeneratedMesh block

  all_GeneratedMesh(1:NumberofGeneratedMesh)%GeneratedMeshId(1)                  = " "
  all_GeneratedMesh(1:NumberofGeneratedMesh)%GeneratedMeshType(1)                = " "
  all_GeneratedMesh(1:NumberofGeneratedMesh)%GeneratedMeshGeometricExtents(1)    = " "
  all_GeneratedMesh(1:NumberofGeneratedMesh)%GeneratedMeshGeometricExtents(2)    = " "
  all_GeneratedMesh(1:NumberofGeneratedMesh)%GeneratedMeshGeometricExtents(3)    = " "
  all_GeneratedMesh(1:NumberofGeneratedMesh)%GeneratedMeshOriginSet(1)           = " "
  all_GeneratedMesh(1:NumberofGeneratedMesh)%GeneratedMeshOriginSet(2)           = " "
  all_GeneratedMesh(1:NumberofGeneratedMesh)%GeneratedMeshOriginSet(3)           = " "
  all_GeneratedMesh(1:NumberofGeneratedMesh)%GeneratedMeshNumberOfElements(1)    = " "
  all_GeneratedMesh(1:NumberofGeneratedMesh)%GeneratedMeshNumberOfElements(2)    = " "
  all_GeneratedMesh(1:NumberofGeneratedMesh)%GeneratedMeshNumberOfElements(3)    = " "

  !! default type for dependnet field block

  all_DependentField(1:NUmberOfDependentFIeld)%DependentFieldId(1)                          = " "
  all_DependentField(1:NUmberOfDependentFIeld)%DependentFieldNumberOfComponents(1)          = " "
  all_DependentField(1:NUmberOfDependentFIeld)%DependentFieldInitialValueOfStateVector(1)   = "0"
  all_DependentField(1:NUmberOfDependentFIeld)%DependentFieldInitialValueOfStateVector(2)   = "0"
  all_DependentField(1:NUmberOfDependentFIeld)%DependentFieldInitialValueOfStateVector(3)   = "0"
  all_DependentField(1:NUmberOfDependentFIeld)%DependentFieldLabel(1)                       = "DEPENDENT_FIELD"
  all_DependentField(1:NUmberOfDependentFIeld)%DependentFieldInitialValueOfStateScalar(1)   = " "


  !! default values for the Region block

  all_Region(1:NumberOfRegion)%RegionId(1)    = " "
  all_Region(1:NumberOfRegion)%RegionLabel(1) = "REGION"

  !!  default values for the Decomposition block

  all_Decomposition(1:NumberOfDecomposition)%CalculateElementFaces(1) = "FALSE"
  all_Decomposition(1:NumberOfDecomposition)%DecompositionId(1)       =  " "

  !!  default values for the Mesh block

  all_Mesh(1:NumberOfMesh)%MeshId(1)   = " "

  !! default value for fields block

  all_Fields(1:NumberOfFields)%FieldsId(1) = " "

  !! default value for 0utput block

  all_Output(1:NUmberOfOutput)%OutputID(1)      = " "
  all_Output(1:NUmberOfOutput)%NodeExport(1)    = "NODE_DATA"
  all_Output(1:NUmberOfOutput)%ElementExport(1) = "ELEMENT_DATA "

  !! default values for source field
  all_SourceField(1:NUmberOfSourceField)%SourceFieldId(1)           = " "
  all_SourceField(1:NUmberOfSourceField)%SourceFieldType(1)         = " "
  all_SourceField(1:NUmberOfSourceField)%SourceFieldComponents(1)   = " "
  all_SourceField(1:NUmberOfSourceField)%SourceFieldComponents(2)   = " "
  all_SourceField(1:NUmberOfSourceField)%SourceFieldComponents(3)   = " "

  !! default value of function

 all_Function(1:NUmberOfFUnction)%FUnctionId(1)                            = " "
 all_Function(1:NUmberOfFUnction)%FUnctionConstants(1)                     = " "
 all_Function(1:NUmberOfFUnction)%FUnctionConstants(2)                     = " "
 all_Function(1:NUmberOfFUnction)%FUnctionConstants(3)                     = " "
 all_Function(1:NUmberOfFUnction)%FUnctionConstants(4)                     = " "
 all_Function(1:NUmberOfFUnction)%FUnctionConstants(5)                     = " "
 all_Function(1:NUmberOfFUnction)%FUnctionConstants(6)                     = " "
 all_Function(1:NUmberOfFUnction)%FUnctionConstants(7)                     = " "
 all_Function(1:NUmberOfFUnction)%FUnctionConstants(8)                     = " "
 all_Function(1:NUmberOfFUnction)%FUnctionConstants(9)                     = " "
 all_Function(1:NUmberOfFUnction)%FUnctionConstants(10)                    = " "

  !! default value for solvers equation block
 !all_SolverEquations(1:NumberOfSolver)%SolverEquationId(1)                  = " "

  !! default value for the fiber field block
 all_FIbreField(1:NumberOfFiberField)%FibreFieldId(1)                     = " "
 all_FIbreField(1:NumberOfFiberField)%FibreFieldParameters(1)             = " "
 all_FIbreField(1:NumberOfFiberField)%FibreFieldParameters(2)             = " "
 all_FIbreField(1:NumberOfFiberField)%FibreFieldParameters(3)             = " "
 all_FibreFIeld(1:NumberOfFiberField)%FiberFieldLabel(1)                  = " "

  !! default value for the fiber field solver

 all_MaterialFIeld(1:NumberOfMaterialField)%MaterialFieldId(1)            = " "
 all_MaterialFIeld(1:NumberOfMaterialField)%MaterialFieldParameters(1)    = " "
 all_MaterialFIeld(1:NumberOfMaterialField)%MaterialFieldParameters(2)    = " "
 all_MaterialFIeld(1:NumberOfMaterialField)%MaterialFieldParameters(3)    = " "
 all_MaterialFIeld(1:NumberOfMaterialField)%MaterialFieldParameters(4)    = " "
 all_MaterialFIeld(1:NumberOfMaterialField)%MaterialFieldParameters(5)    = " "
 all_MaterialFIeld(1:NumberOfMaterialField)%MaterialFieldParameters(6)    = " "
 all_MaterialFIeld(1:NumberOfMaterialField)%MaterialFieldParameters(7)    = " "
 all_MaterialFIeld(1:NumberOfMaterialField)%MaterialFieldParameters(8)    = " "
 all_MaterialFIeld(1:NumberOfMaterialField)%MaterialFieldParameters(9)    = " "
 all_MaterialFIeld(1:NumberOfMaterialField)%MaterialFieldParameters(10)   = " "
 all_MaterialFIeld(1:NumberOfMaterialField)%MaterialFieldLabel(1)         = "MATERIAL"

