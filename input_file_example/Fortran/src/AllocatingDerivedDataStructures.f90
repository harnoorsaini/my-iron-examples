
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!! THE FOLLOWING HEADER FILE PERFOMS THE FOLLOWINGS TASKS. ... !!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 1) PARSE THROUGH THE INPUT FILE TO COUNT AND STORE BASIS, REGIONS ETC. DEFINED... !!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ...BY USER IN THE INPUT FILE.  ...!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 2) BASED ON THE COUNTS IT ALLOCATES THE DERIVED DATA STRUCTURE. !!!!!!!!!!!!!!!!!!!!!


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!! FOLLOWING VARIABLES STORE  NUMBER OF TIMES BASIS , PARTS , REGIONS ETC.  .....!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!! .... ARE DEFINED BY THE USER. !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  NumberOfEquationsSet= 0
  NumberOfProblem = 0
  NumberOfBasis = 0
  NumberOfSolver = 0
  NumberOfMesh = 0
  NumberOfBoundaryCondition = 0
  NumberOfBasis = 0
  NumberOfMaterialField = 0
  NumberOfEquation = 0
  NumberOfGeometricField = 0
  NumberOfControlLoop = 0
  NumberOfDependentField = 0
  NumberOfRegion = 0
  NumberOfCoordinateSystem = 0
  NumberOfFiberField       = 0
  NumberOfSourceField    = 0
  NumberOfOutput   = 0
  NumberOfDecomposition = 0
  NumberOfGeneratedMesh = 0
  NumberOfFields = 0
  NumberOfFunction = 0


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!! PARSE INPUT FILE TO STORE INFORMATION. !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  OPEN(12,file=InputFile ,status="old")
  READ(12,'(A)') rdline

  DO WHILE (TRIM(rdline).NE."STOP_PARSING")

       READ(12,'(A)') rdline
       CALL searching(rdline,"START_EQUATIONS_SET",NumberOfEquationsSet)
       CALL searching(rdline,"START_PROBLEM",NumberOfProblem)
       CALL searching(rdline,"START_BOUNDARY_CONDITIONS",NumberOfBoundaryCondition)
       CALL searching(rdline,"START_BASIS",NumberOfBasis)
       CALL searching(rdline,"START_MESH",NumberOfMesh)
       CALL searching(rdline,"START_GENERATED_MESH",NumberOfGeneratedMesh)
       CALL searching(rdline,"START_MATERIAL_FIELD",NumberOfMaterialField)
       CALL searching(rdline,"START_REGION",NumberOfRegion)
       CALL searching(rdline,"START_GEOMETRIC_FIELD",NumberOfGeometricField)
       CALL searching(rdline,"START_DEPENDENT_FIELD",NumberOfDependentField)
       CALL searching(rdline,"START_SOLVER_SETTINGS",NumberOfSolver)
       CALL searching(rdline,"START_CONTROL_LOOP",NumberOfControlLoop)
       CALL searching(rdline,"START_FIBER_FIELD",NumberOfFiberField)
       CALL searching(rdline,"START_SOURCE_FIELD",NumberOfSourceField)
       CALL searching(rdline,"START_OUTPUT",NumberOfOutput)
       CALL searching(rdline,"START_COORDINATE_SYSTEM",NumberOfCoordinateSystem)
       CALL searching(rdline,"START_DECOMPOSITION",NumberOfDecomposition)
       CALL searching(rdline,"START_FIELDS",NumberOfFields)
       CALL searching(rdline,"START_FUNCTION",NumberOfFunction)

  END DO

  !!!NumberOfGeneratedMesh = NumberOfMesh    	     !! PLEASE IGNORE IT, I WILL FIX IT AFTER OUR MEETING
  !!!!!!!!!!!!!!!!!		ALLOCATE DATA STRUCTURES BASED ON THE PARAMETERS INITIALIZED ABOVE.  !!!!!!!!!!!!!!!!!

  IF (NumberOfBasis .GT. 0 )   ALLOCATE(all_Basis(NumberOfBasis))

  IF (NumberOfBoundaryCondition .GT. 0 )   ALLOCATE(all_BoundaryConditions(NumberOfBoundaryCondition))

  IF (NumberOfCoordinateSystem .GT. 0 )   ALLOCATE(all_CoordinateSystem(NumberOfCoordinateSystem))

  IF (NumberOfWorldCoordinateSystem .GT. 0 )    &
    & ALLOCATE(all_WorldCoordinateSystem%WorldCoordinateSystem(NumberOfWorldCoordinateSystem))

  IF (NumberOfMesh .GT. 0 )   ALLOCATE(all_Mesh(NumberOfMesh))

  IF (NumberOfDecomposition .GT. 0 )   ALLOCATE(all_Decomposition(NumberOfDecomposition))

  IF (NumberOfEquationsSet .GT. 0 )   ALLOCATE(all_Equations%Equations(NumberOfEquationsSet))

  IF (NumberOfEquationsSet .GT. 0 )   ALLOCATE(all_EquationsSet(NumberOfEquationsSet))

  IF (NumberOfGeometricField .GT. 0 )   ALLOCATE(all_GeometricField(NumberOfGeometricField))

  IF (NumberOfFiberField .GT. 0 )   ALLOCATE(all_FibreField(NumberOfFiberField))

  IF (NumberOfMaterialField .GT. 0 )   ALLOCATE(all_MaterialField(NumberOfMaterialField))

  IF (NumberOfEquationsSet .GT. 0 )   ALLOCATE(all_EquationsSetField(NumberOfEquationsSet))

  IF (NumberOfFields .GT. 0 )   ALLOCATE(all_Fields(NumberOfFields))

  IF (NumberOfProblem .GT. 0 )   ALLOCATE(all_Problem(NumberOfProblem))

  IF (NumberOfRegion .GT. 0 )   ALLOCATE(all_Region(NumberOfRegion))

  IF (NumberOfWorldRegion .GT. 0 )   ALLOCATE(all_WorldRegion%WorldRegion(NumberOfWorldRegion))

  IF (NumberOfSolver .GT. 0 )   ALLOCATE(all_Solver(NumberOfSolver))

  IF (NumberOfSolver .GT. 0 )   ALLOCATE(all_LinearSolver(NumberOfSolver))

  IF (NumberOfSolver .GT. 0 )   ALLOCATE(all_NonLinearSolver(NumberOfSolver))

  IF (NumberOfSolver .GT. 0 )   ALLOCATE(all_SolverEquations(NumberOfSolver))

  IF (NumberOfControlLoop .GT. 0 )   ALLOCATE(all_ControlLoop(NumberOfControlLoop))

  IF (NumberOfGeneratedMesh .GT. 0 )   ALLOCATE(all_GeneratedMesh(NumberOfGeneratedMesh))

  IF (NumberOfDependentField .GT. 0 )   ALLOCATE(all_DependentField(NumberOfDependentField))

  IF (NumberOfSourceField .GT. 0 )   ALLOCATE(all_SourceField(NumberOfSourceField))

  IF (NumberOfOutput .GT. 0 )   ALLOCATE(all_Output(NumberOfOutput))

  IF (NumberOfFunction .GT. 0 )   ALLOCATE(all_Function(NumberOfFunction))

