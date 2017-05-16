
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!! _______________________ HEADER COMMENTS ____________________!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!! THE FOLLOWING HEADER FILE INITILIZES THE POINTER TYPE DERIVED DATA STRCUTRES....!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!! .... WITH DERIVED TYPES. !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!________________________________________________________!!!!!!!!!!!!!!!!!!!!!!!!



  !!!!!!!!!!!!!!!!!!!!!!!!! DEFINING DERIVED TYPES AND ALLOCATING THIER RESPECTIVE DATA STRUCTURES !!!!!!!!!!!!!!!!!!!

  TYPE BasisStructure

    TYPE(CMFE_BasisType) 	                       :: Basis
    TYPE(CHARACTER(LEN=100))                           :: BasisId(1),BasisInterpolationType(1), BasisNumberOfGaussPoints(3)

  END TYPE BasisStructure

  TYPE BoundaryConditionStructure

    TYPE(CMFE_BoundaryConditionsType)                  :: BoundaryConditions


  END TYPE  BoundaryConditionStructure

  TYPE CoordinateSystemStructure

    TYPE(CMFE_CoordinateSystemType)    , ALLOCATABLE :: WorldCoordinateSystem(:)
    TYPE(CHARACTER(LEN=100))                         :: CoordinateSystemId(1),CoordinateSystemType(1),CoordinateSystemDimension(1)
    TYPE(CMFE_CoordinateSystemType)                  :: CoordinateSystem

  END TYPE CoordinateSystemStructure

  TYPE MeshStructure

    TYPE(CMFE_MeshType)                                :: Mesh
    TYPE(CHARACTER(LEN=100))                           :: MeshId(1)

  END TYPE  MeshStructure

  TYPE DecompositionStructure

    TYPE(CMFE_DecompositionType)                      :: Decomposition
    TYPE(CHARACTER(LEN=100))                          :: CalculateElementFaces(1),DecompositionId(1)

  END TYPE  DecompositionStructure

  TYPE equations_array

    TYPE(CMFE_EquationsType)            , ALLOCATABLE :: Equations(:)

  END TYPE equations_array

  TYPE EquationsSetStructure

    TYPE(CMFE_EquationsSetType)                        :: EquationsSet
    TYPE(CHARACTER(LEN=100))                           :: EquationSetId(1),EquationSetClass(1),EquationSetType(1), &
                                                          & EquationSetSubType(1),EquationSetOutputTypeSet(1)

  END TYPE EquationsSetStructure

  TYPE FieldTypeStructure

    TYPE(CMFE_FieldType)                               :: GeometricField,MaterialField,EquationsSetField,FibreField
    TYPE(CHARACTER(LEN=100))                           :: GeometricFieldId(1),GeometricFieldLabel(1), EquationsSetFieldId(1)
    TYPE(CHARACTER(LEN=100))                           :: FibreFieldId(1),FibreFieldParameters(3),FiberFieldLabel(1)
    TYPE(CHARACTER(LEN=100))                           :: MaterialFieldId(1),MaterialFieldParameters(10),MaterialFieldLabel(1)

  END TYPE  FieldTypeStructure

  TYPE FieldsStructure

    TYPE(CMFE_FieldsType)                              ::  Fields
    TYPE(CHARACTER(LEN=100))                           ::  FieldsId(1)

  END TYPE  FieldsStructure

  TYPE ProblemStructure

    TYPE(CMFE_ProblemType)                             :: Problem
    CHARACTER(LEN=100)                                 :: ProblemId(1),ProblemClass(1),ProblemType(1),ProblemSubType(1)

  END TYPE  ProblemStructure

  TYPE RegionStructure

    TYPE(CMFE_RegionType)                              :: Region
    CHARACTER(LEN=100)                                 :: RegionId(1),RegionLabel(1)
    TYPE(CMFE_RegionType)                , ALLOCATABLE :: WorldRegion(:)

  END TYPE RegionStructure

  TYPE SolverStructure

    TYPE(CMFE_SolverType)                             :: Solver,LinearSolver,NonLinearSolver
    TYPE(CHARACTER(LEN=100))                          :: SolverId(1) , LinearSolverType(1), DirectSolverType(1)
    TYPE(CHARACTER(LEN=100))                          :: IterativeSolverType(1),IterativeSolverRelativeTolerance(1)
    TYPE(CHARACTER(LEN=100))                          :: IterativeSolverAbsoluteTolerance(1),IterativeSolverMaximumIterations(1)
    TYPE(CHARACTER(LEN=100))                          :: IterativePreconditioner(1),NonlinearSolverType(1)
    TYPE(CHARACTER(LEN=100))                          :: NewtonMaximumIterations(1),NewtonRelativeTolerance(1)
    TYPE(CHARACTER(LEN=100))                          :: NewtonAbsoluteTolerance(1),OutputTypeSet(1)
    TYPE(CHARACTER(LEN=100))                          :: JacobianType(1) , LibraryType(1)

  END TYPE SolverStructure

  TYPE SolverEquationsStructure

    TYPE(CMFE_SolverEquationsType)                     :: SolverEquations

  END TYPE SolverEquationsStructure

  TYPE ControlLoopStructure

    TYPE(CMFE_ControlLoopType)                         :: ControlLoop
    TYPE(CHARACTER(LEN=100))                           :: ControlLoopId(1),ControlLoopType(1),ControlLoopTimeIncrement(3), &
                                                          & ControlLoopLoadIncrement(3)

  END TYPE ControlLoopStructure

  TYPE GeneratedMeshStructure

    TYPE(CMFE_GeneratedMeshType)                       :: GeneratedMesh
    TYPE(CHARACTER(LEN=100))                           :: GeneratedMeshID(1),GeneratedMeshType(1)
    TYPE(CHARACTER(LEN=100))                           :: GeneratedMeshGeometricExtents(3),GeneratedMeshOriginSet(3)
    TYPE(CHARACTER(LEN=100))                           :: GeneratedMeshNumberOfElements(3)

  END TYPE GeneratedMeshStructure

  TYPE DependentFieldStructure

    TYPE(CMFE_FieldType)                              :: DependentField
    TYPE(CHARACTER(LEN=100))                          :: DependentFieldId(1),DependentFieldStateVariable(1)
    TYPE(CHARACTER(LEN=100))                          :: DependentFieldNumberOfComponents(1)
    TYPE(CHARACTER(LEN=100))                          :: DependentFieldInitialValueOfStateVector(3), DependentFieldLabel(1)
    TYPE(CHARACTER(LEN=100))                          :: DependentFieldInitialValueOfStateScalar(1)


  END TYPE DependentFieldStructure


  TYPE SourceFieldStructure

    TYPE(CMFE_FieldType)                               :: SourceField
    TYPE(CHARACTER(LEN=100))                           :: SourceFieldId(1) , SourceFieldType(1), SourceFieldComponents(3)


  END TYPE SourceFieldStructure


  TYPE OutputStructure

    TYPE(CHARACTER(LEN=100))                           :: OutputID(1),NodeExport(1),ElementExport(1)

  END TYPE OutputStructure


  TYPE FunctionStructure

    TYPE(CHARACTER(LEN=100))                           :: FunctionID(1),FunctionConstants(10)

  END TYPE FunctionStructure

 !!!!!! INITIALIZING POINTER THAT POINT AT THE DATA STRUCTURES OF DERIVED TYPES !!!!!!!!!!!!!!!!

 TYPE(BasisStructure)             , ALLOCATABLE    :: all_Basis(:)
 TYPE(BoundaryConditionStructure) , ALLOCATABLE    :: all_BoundaryConditions(:)
 TYPE(CoordinateSystemStructure)  , ALLOCATABLE    :: all_CoordinateSystem(:)
 TYPE(CoordinateSystemStructure)                   :: all_WorldCoordinateSystem
 TYPE(MeshStructure)              , ALLOCATABLE    :: all_Mesh(:)
 TYPE(DecompositionStructure) 	  , ALLOCATABLE    :: all_Decomposition(:)
 TYPE(equations_array) 	                           :: all_Equations
 TYPE(EquationsSetStructure)      , ALLOCATABLE    :: all_EquationsSet(:)
 TYPE(FieldTypeStructure)         , ALLOCATABLE    :: all_GeometricField(:)
 TYPE(FieldTypeStructure)         , ALLOCATABLE    :: all_FibreField(:)
 TYPE(FieldTypeStructure)         , ALLOCATABLE    :: all_MaterialField(:)
 TYPE(FieldTypeStructure)         , ALLOCATABLE    :: all_EquationsSetField(:)
 TYPE(FieldsStructure)            , ALLOCATABLE    :: all_Fields(:)
 TYPE(ProblemStructure)           , ALLOCATABLE    :: all_Problem(:)
 TYPE(RegionStructure)            , ALLOCATABLE    :: all_Region(:)
 TYPE(RegionStructure)		                   :: all_WorldRegion
 TYPE(SolverStructure)		  , ALLOCATABLE    :: all_Solver(:)
 TYPE(SolverStructure)		  , ALLOCATABLE    :: all_LinearSolver(:)
 TYPE(SolverStructure)		  , ALLOCATABLE    :: all_NonLinearSolver(:)
 TYPE(SolverEquationsStructure)   , ALLOCATABLE    :: all_SolverEquations(:)
 TYPE(ControlLoopStructure)       , ALLOCATABLE    :: all_ControlLoop(:)
 TYPE(GeneratedMeshStructure)     , ALLOCATABLE    :: all_GeneratedMesh(:)
 TYPE(DependentFieldStructure)    , ALLOCATABLE    :: all_DependentField(:)
 TYPE(SourceFieldStructure)       , ALLOCATABLE    :: all_SourceField(:)
 TYPE(OUtputStructure)            , ALLOCATABLE    :: all_Output(:)
 TYPE(FunctionStructure)          , ALLOCATABLE    :: all_Function(:)

