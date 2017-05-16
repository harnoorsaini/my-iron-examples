

 !!!!!!!!!!!!!!!!!!!!________THE FOLLOWING HEADER FILE ALLOCATES THE SIZE OF "USER NUMBER" DATA STRUCTURES BASED ON THE VARIABLES INITIALIED IN _____!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!! ....... HEADER FILE ALLOCATING DERIVED DATA STRUCTURES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!________________________________________________________________________________________________________________!!!!!!!!!!!!!!!!

  ALLOCATE (CoordinateSystemUserNumber(NumberOfCoordinateSystem))
  ALLOCATE (RegionUserNumber(NumberOfRegion))
  ALLOCATE (BasisUserNumber(NumberOfBasis))
  ALLOCATE (GeneratedMeshUserNumber(NumberOfGeneratedMesh))
  ALLOCATE (MeshUserNumber(NumberOfMesh))
  ALLOCATE (DecompositionUserNumber(NumberOfDecomposition))
  ALLOCATE (FieldGeometryUserNumber(NumberOfGeometricField))
  ALLOCATE (FieldFibreUserNumber(NumberOfFiberField))
  ALLOCATE (FieldMaterialUserNumber(NumberOfMaterialField))
  ALLOCATE (FieldDependentUserNumber(NumberOfDependentField))
  ALLOCATE (EquationSetUserNumber(NumberOfEquationsSet))
  ALLOCATE (EquationsSetFieldUserNumber(NumberOfEquationsSet))
  ALLOCATE (ProblemUserNumber(NumberOfProblem))
  ALLOCATE (FieldUserNumber(NumberOfProblem))


!!!!!!!!!!!!!!!!!! Assigning tags. NOTE:  Each tag should be distinct !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 UserNumberLabel = 1

  DO i = 1, NumberOfCoordinateSystem

        CoordinateSystemUserNumber(i)  = UserNumberLabel
        UserNumberLabel = UserNumberLabel + 1
  END DO

 DO i = 1, NumberOfRegion

        RegionUserNumber(i)  = UserNumberLabel
        UserNumberLabel = UserNumberLabel + 1
 END DO


 DO i = 1, NumberOfBasis

        BasisUserNumber(i)  = UserNumberLabel
        UserNumberLabel = UserNumberLabel + 1
 END DO


 DO i = 1,NumberOfGeneratedMesh

        GeneratedMeshUserNumber(i)  = UserNumberLabel
        UserNumberLabel = UserNumberLabel + 1
 END DO

 DO i = 1, NumberOfMesh

       MeshUserNumber(i) = UserNumberLabel
       UserNumberLabel = UserNumberLabel + 1
 END DO

 DO i = 1, NumberOfDecomposition

       DecompositionUserNumber(i) = UserNumberLabel
       UserNumberLabel = UserNumberLabel + 1
 END DO

 DO i = 1, NumberOfGeometricField

       FieldGeometryUserNumber(i) = UserNumberLabel
       UserNumberLabel = UserNumberLabel + 1
 END DO

 DO i = 1, NumberOfFiberField

       FieldFibreUserNumber(i) = UserNumberLabel
       UserNumberLabel = UserNumberLabel + 1
 END DO

 DO i = 1, NumberOfMaterialField

       FieldMaterialUserNumber(i) = UserNumberLabel
       UserNumberLabel = UserNumberLabel + 1
 END DO

 DO i = 1, NumberOfDependentField

       FieldDependentUserNumber(i) = UserNumberLabel
       UserNumberLabel = UserNumberLabel + 1
 END DO

 DO i = 1, NumberOfEquationsSet

       EquationSetUserNumber(i) = UserNumberLabel
       UserNumberLabel = UserNumberLabel + 1

 END DO

 DO i = 1, NumberOfEquationsSet

       EquationsSetFieldUserNumber(i) = UserNumberLabel
       UserNumberLabel = UserNumberLabel + 1
 END DO

 DO i = 1, NumberOfProblem

       ProblemUserNumber(i) = UserNumberLabel
       UserNumberLabel = UserNumberLabel + 1
 END DO

 DO i = 1, NumberOfSourceField

        FieldUserNumber(i)  = UserNumberLabel
        UserNumberLabel = UserNumberLabel + 1
 END DO

