!> \file
!> \author: Waleed Mirza
!> \brief This is a generic OpenCMISS-iron example program to solve various different types of equations using OpenCMISS CALLs.
!>
!> \section LICENSE
!>
!> Version: MPL 1.1/GPL 2.0/LGPL 2.1
!>
!> The contents of this file are subject to the Mozilla Public License
!> Version 1.1 (the "License"); you may not use this file except in
!> compliance with the License. You may obtain a copy of the License at
!> http://www.mozilla.org/MPL/
!>
!> Software distributed under the License is distributed on an "AS IS"
!> basis, WITHOUT WARRANTY OF ANY KIND, either express or implied. See the
!> License for the specific language governing rights and limitations
!> under the License.
!>
!> The Original Code is OpenCMISS
!>
!> The Initial Developer of the Original Code is University of Auckland,
!> Auckland, New Zealand and University of Oxford, Oxford, United
!> \file
!> \author: Waleed Mirza
!> \brief This is a generic OpenCMISS-iron example program to solve various different types of equations using OpenCMISS CALLs.
!>
!> \section LICENSE
!>
!> Version: MPL 1.1/GPL 2.0/LGPL 2.1
!>
!> The contents of this file are subject to the Mozilla Public License
!> Version 1.1 (the "License"); you may not use this file except in
!> compliance with the License. You may obtain a copy of the License at
!> http://www.mozilla.org/MPL/
!>
!> Software distributed under the License is distributed on an "AS IS"
!> basis, WITHOUT WARRANTY OF ANY KIND, either express or implied. See the
!> License for the specific language governing rights and limitations
!> under the License.
!>
!> The Original Code is OpenCMISS
!>
!> The Initial Developer of the Original Code is University of Auckland,
!> Auckland, New Zealand and University of Oxford, Oxford, United
!> KingDOm. Portions created by the University of Auckland and University
!> of Oxford are Copyright (C) 2007 by the University of Auckland and
!> the University of Oxford. All Rights Reserved.
!>
!> Contributor(s): Waleed Mirza, Andreas Hessenthaler, Thomas Heidlauf, Thomas Klotz
!>
!> Alternatively, the contents of this file may be used under the terms of
!> either the GNU General Public License Version 2 or later (the "GPL"), or
!> the GNU Lesser General Public License Version 2.1 or later (the "LGPL"),
!> in which case the provisions of the GPL or the LGPL are applicable instead
!> of those above. IF you wish to allow use of your version of this file only
!> under the terms of either the GPL or the LGPL, and not to allow others to
!> use your version of this file under the terms of the MPL, indicate your
!> decision by deleting the provisions above and replace them with the notice
!> and other provisions required by the GPL or the LGPL. IF you DO not delete
!> the provisions above, a recipient may use your version of this file under
!> the terms of any one of the MPL, the GPL or the LGPL.
!>


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! THE FOLLOWING HEADER FILE INCLUDES THE DATA STRUCTURES  SUBROUTINES AND FUNCTIONS EMPLOYED IN PARSING ALGORITHM.  !!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! FOR MORE INFO READ HEADER COMMENTS OF THE FILE. !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
INCLUDE "parsing_module.f90"


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! MAIN PROGRAM STARTS FROM HERE. !!!!!!!!!!!!!!!!!!
PROGRAM GENERICEXAMPLE

  USE OpenCMISS
  USE OpenCMISS_Iron
  USE Parsing

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! LOGICAL PREPROCESSOR  DIRECTIVES. !!!!!!!!!!!!!!!
#ifndef NOMPIMOD
  USE MPI
#endif

#ifdef WIN32
  USE IFQWIN
#endif

  IMPLICIT NONE

#ifdef NOMPIMOD
#INCLUDE "mpif.h"
#endif


  !!!!!!!!!!!!!!!!!!!!VARIABLES DEFINING SIZES OF DERIVED DATA STRUCTURES. !!!!!!!!!!!!!!!

  INTEGER(CMISSINTG)           :: NumberOfBasis
  INTEGER(CMISSINTG)           :: NumberOfBoundaryCondition
  INTEGER(CMISSINTG),PARAMETER :: NumberOfWorldCoordinateSystem = 1 ! ALways equal to 1.
  INTEGER(CMISSINTG),PARAMETER :: NumberOfWorldRegion = 1           ! Always equal to 1.
  INTEGER(CMISSINTG)           :: NumberOfCoordinateSystem
  INTEGER(CMISSINTG)           :: NumberOfMesh
  INTEGER(CMISSINTG)           :: NumberOfDecomposition
  INTEGER(CMISSINTG)           :: NumberOfEquation
  INTEGER(CMISSINTG)           :: NumberOfEquationsSet
  INTEGER(CMISSINTG)           :: NumberOfGeometricField
  INTEGER(CMISSINTG)           :: NumberOfFiberField
  INTEGER(CMISSINTG)           :: NumberOfMaterialField
  INTEGER(CMISSINTG)           :: NumberOfDependentField
  INTEGER(CMISSINTG)           :: NumberOfEquationSetField
  INTEGER(CMISSINTG)           :: NumberOfSourceField
  INTEGER(CMISSINTG)           :: NumberOfProblem
  INTEGER(CMISSINTG)           :: NumberOfRegion
  INTEGER(CMISSINTG)           :: NumberOfSolver
  INTEGER(CMISSINTG)           :: NumberOfLinearSolver
  INTEGER(CMISSINTG)           :: NumberOfSolverEquations
  INTEGER(CMISSINTG)           :: NumberOfControlLoop
  INTEGER(CMISSINTG)           :: NumberOfGeneratedMesh
  INTEGER(CMISSINTG)           :: NumberOfOutput
  INTEGER(CMISSINTG)           :: NumberOfFields
  INTEGER(CMISSINTG)           :: NumberOfFunction

  !!  These store the number of boundary conditions defined by user in the input file.
  INTEGER(CMISSINTG)                :: NumberOfDirichelet,NumberOfTractionNeumann,NumberOfPressureNeumann
  !!!!!!!!!!!!!!!!!!!!        COUNTERS USED IN FORTRAN EXAMPLE FILE.                    !!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!   FOR INSTANCE "SolverIdx" IS USED IN THE BLOCK ........     !!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!  ......WHERE SOLVER IS DEFINED.                                  !!! !!!!!!!!!!!!

  INTEGER(CMISSINTG)          :: CoordinateSystemIdx, RegionIdx,BasisIdx,DecompositionIdx,BoundaryConditionIdx,SourceFieldIdx
  INTEGER(CMISSINTG)          :: ProblemIdx,SolverIdx,EquationSetIdx,MaterialFieldIdx,DependentFieldIdx,FunctionIdx
  INTEGER(CMISSINTG)          :: FiberFieldIdx,GeometricFieldIdx,GeneratedMeshIdx,ControlLoopIdx,FieldIdx,MeshIdx,SecondBasisIdx
  !! Counters used in different sections of FortranExample.f90.
  INTEGER(CMISSINTG)          :: j,k ,h, UserNumberLabel, ComponentIdx,MaterialParametersIdx
  INTEGER(CMISSINTG)          :: ProblemId, BasisId , BoundaryId, OutputId, RegionId,EquationId,FieldsId, FunctionId
  INTEGER(CMISSINTG)          :: ControlLoopId, SolverSettingsId , DependentFieldId, CoordinateSystemId, DecompositionId
  INTEGER(CMISSINTG)          :: MaterialFieldId, FibreFieldId , MeshId,GeneratedMeshId, SourceFieldId,GeometricFieldId
  !! Parameter to be used in Boundary Condition bit.
  CHARACTER(LEN = 100)        :: Constraint,FunctionName

  !!!!!!!!!!!!!!!! FOLLOWING DATA STRUCTURE STORE THE ID´s ASSINGED TO REGIONS, BASIS etc. !!!!!!!!

  INTEGER(CMISSINTG), ALLOCATABLE :: CoordinateSystemUserNumber(:)
  INTEGER(CMISSINTG), ALLOCATABLE :: RegionUserNumber(:)
  INTEGER(CMISSINTG), ALLOCATABLE :: BasisUserNumber(:)
  INTEGER(CMISSINTG), ALLOCATABLE :: GeneratedMeshUserNumber(:)
  INTEGER(CMISSINTG), ALLOCATABLE :: MeshUserNumber(:)
  INTEGER(CMISSINTG), ALLOCATABLE :: DecompositionUserNumber(:)
  INTEGER(CMISSINTG), ALLOCATABLE :: FieldGeometryUserNumber(:)
  INTEGER(CMISSINTG), ALLOCATABLE :: FieldFibreUserNumber(:)
  INTEGER(CMISSINTG), ALLOCATABLE :: FieldMaterialUserNumber(:)
  INTEGER(CMISSINTG), ALLOCATABLE :: FieldDependentUserNumber(:)
  INTEGER(CMISSINTG), ALLOCATABLE :: EquationSetUserNumber(:)
  INTEGER(CMISSINTG), ALLOCATABLE :: EquationsSetFieldUserNumber(:)
  INTEGER(CMISSINTG), ALLOCATABLE :: ProblemUserNumber(:)
  INTEGER(CMISSINTG), ALLOCATABLE :: FieldUserNumber(:)
  INTEGER(CMISSINTG), ALLOCATABLE :: BasisIdxArray(:,:)


  !!!!!!!!!!!!!!! FOLLOWING ARE THE PROGRAM VARIABLES  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  INTEGER(CMISSINTG)              :: Err
  REAL(CMISSRP)                   :: SecondDimension,FirstDImension,ThirdDImension                              !! These variables stores dimensions of the DOmain.
  INTEGER(CMISSINTG)              :: NumberGlobalXElements,NumberGlobalYElements,NumberGlobalZElements          !! These variables stores  number of elements in each dimension.
  INTEGER(CMISSINTG)              :: NumberOfArguments,ArgumentThirdDImension,ArgStatus                         !! These variables are used in reading input arguments "path + file name"" ..
  CHARACTER(LEN=300)              :: CommandArgument,InputFile                                                  !! ... from the command line. 
  INTEGER(CMISSINTG)              :: InterpolationType,NumberOfGaussXi(3),NumberOfComponents,DependentFieldComponent !! These variables store information to be used in "BASIS BLOCK"  defined by ..
  INTEGER(CMISSINTG)              :: ScalingType                                                                !! As of now, i have no idea what this is for , i will look into it later.
  INTEGER(CMISSINTG)              :: EquationsSetIndex                                                          !! As of now, i have no idea what this is for , i will look into it later.
  INTEGER(CMISSINTG)              :: NumberOfComputationalNodes,NumberOfDomains,ComputationalNodeNumber         !! THese variables are used in "DECOMPOSITION BLOCK".
  INTEGER(CMISSINTG)              :: NodeNumber,NodeDOmain,nodeIdx                                              !! Same as the last comment
  INTEGER(CMISSINTG),ALLOCATABLE  :: SurfaceNodes(:)                                                            !! Possible surface where Dirichelet  or Natural BCs are imposed.

  INTEGER(CMISSIntg)              ::  LeftNormalXi

                                 !! Normal directive of possible surfaces where Dirichelet  or Natural BCs are imposed.
  INTEGER(CMISSINTG),ALLOCATABLE  :: PressureComponentExist(:)
  REAL(CMISSRP)                   :: x,y,z                                                                      !! spatial coordinates
  REAL(CMISSRP)                   :: BOundaryConditionValue

  !!!!!!!!!!!!!!!!!!!!!!!!!!! FOLLOWING HEADER FILE INITIALIZE  THE DATA STRUCTURES WITH DERIVED TYPES.!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! FOR MORE INFO READ HEADER COMMENTS OF THE FILE. !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  INCLUDE "DerivedTypes.f90"

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! LOGICAL PREPROCESSOR  DIRECTIVES. !!!!!!!!!!!!!!!
#ifdef WIN32
  !Quickwin type
  LOGICAL :: QUICKWIN_STATUS=.FALSE.
  TYPE(WINDOWCONFIG) :: QUICKWIN_WINDOW_CONFIG
#endif


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! THE FOLLOWING BLOCK READS THE PATH AND NAME OF THE INPUT FILE.   !!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! IT MAKES USE OF THE INTRINSIC FUNCTIONS. : ...................   !!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ... 1- COMMAND_ARGUMENT_COUNT(). 2- GET_COMMAND_ARGUMENT().      !!!!!!!!!!!!!!!!!!!!!!!!

  NumberOfArguments = COMMAND_ARGUMENT_COUNT()

  IF(NumberOfArguments .NE. 1) THEN

    CALL HANDLE_ERROR("Please provide only the input file")                          !! Throws warning IF the path is wrong.

  ELSE

    CALL GET_COMMAND_ARGUMENT(1,CommandArgument,ArgumentThirdDImension,ArgStatus)
    InputFile = TRIM(CommandArgument)

  END IF



  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  FOLLLOWING HEADER FILE INITIALIZE THE SIZE OF DERIVED DATA STRUCTURES. !!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! FOR MORE INFO READ HEADER COMMENTS OF THE FILE. !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  INCLUDE "AllocatingDerivedDataStructures.f90"

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! FOLLOWING HEADER FILE ALLOCATES  "USERNUMBER"  DATA STRUCTURE SIZE           !!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! FOR MORE INFO READ HEADER COMMENTS OF THE FILE. !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  INCLUDE "UserNumberDataStructures.f90"

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! THE FOLLOWING FILE INITILIZES DEFAULT VALUES OF THE INPUT ARGUMENTS !!!!!!!!!!!!!!!!!!!!!!!
  INCLUDE "DefaultArgument.f90"
    
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        INCLUDES I/0 STATEMENT THAT PARSE THROUGH THE INPUT FILE        !!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! FOR MORE INFO READ HEADER COMMENTS OF THE FILE. !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  INCLUDE "parsing_algorithm.f90"

#ifdef WIN32
  !Initialise QuickWin
  QUICKWIN_WINDOW_CONFIG%TITLE="General Output" !Window title
  QUICKWIN_WINDOW_CONFIG%NUMTEXTROWS=-1 !Max possible number of rows
  QUICKWIN_WINDOW_CONFIG%MODE=QWIN$SCROLLDOWN
  !Set the window parameters
  QUICKWIN_STATUS=SETWINDOWCONFIG(QUICKWIN_WINDOW_CONFIG)
  !If attempt fails set with system estimated values
  IF(.NOT.QUICKWIN_STATUS) QUICKWIN_STATUS=SETWINDOWCONFIG(QUICKWIN_WINDOW_CONFIG)
#endif

  !Intialise cmiss
  CALL cmfe_Initialise(all_WorldCoordinateSystem%WorldCoordinateSystem(1),all_WorldRegion%WorldRegion(1),Err)

  CALL cmfe_ErrorHandlingModeSet(CMFE_ERRORS_TRAP_ERROR,Err)
   !Get the number of computational nodes and this computational node number

  CALL cmfe_ComputationalNumberOfNodesGet(NumberOfComputationalNodes,Err)

  CALL cmfe_ComputationalNodeNumberGet(ComputationalNodeNumber,Err)

  NumberOfDomains=NumberOfComputationalNodes


  !!!!!!!!!!!!! INTRODUCING  LOOPS FOR ASSIGNING COORDINATE SYSTEM TO THIER RESPECTIVE REGIONS !!!!!!!!!!

  DO  CoordinateSystemIdx = 1,NumberOfCoordinateSystem
    !Create a 3D rectangular cartesian coordinate system

    CALL cmfe_CoordinateSystem_Initialise(all_CoordinateSystem(CoordinateSystemIdx)%CoordinateSystem,Err)

    CALL cmfe_CoordinateSystem_CreateStart(CoordinateSystemUserNumber(CoordinateSystemIdx), &
      & all_CoordinateSystem(CoordinateSystemIdx)%CoordinateSystem,Err)

    CALL cmfe_CoordinateSystem_TypeSet(CoordinateSystemUserNumber(CoordinateSystemIdx), &
      & MATCH_COORDINATE_SYSTEM(all_CoordinateSystem(CoordinateSystemIdx)%CoordinateSystemType(1)), err)

    CALL cmfe_CoordinateSystem_DimensionSet(all_CoordinateSystem(CoordinateSystemIdx)%CoordinateSystem, &
      & MATCH_COORDINATE_SYSTEM(all_CoordinateSystem(CoordinateSystemIdx)%CoordinateSystemDimension(1)),Err)

    CALL cmfe_CoordinateSystem_CreateFinish(all_CoordinateSystem(CoordinateSystemIdx)%CoordinateSystem,Err)

  END DO ! CoordinateSystemIdx

  
  DO RegionIdx = 1, NumberOfRegion

    !atch IDs for regions and coordinate system objects
    CALL MATCH_IDs(all_Region(RegionIdx)%RegionID(1),all_CoordinateSystem(:)%CoordinateSystemID(1),CoordinateSystemIdx, &
      & "all_Region(:)", "all_CoordinateSystem(:)" )

    !Create a region and assign the coordinate system to the region
    CALL cmfe_Region_Initialise(all_Region(RegionIdx)%Region,Err)

    CALL cmfe_Region_CreateStart(RegionUserNumber(RegionIdx),all_WorldRegion%WorldRegion(RegionIdx), &
      & all_Region(RegionIdx)%Region,Err)

    CALL cmfe_Region_LabelSet(all_Region(RegionIdx)%Region,all_Region(RegionIdx)%RegionLabel(1),Err)

    CALL cmfe_Region_CoordinateSystemSet(all_Region(RegionIdx)%Region, &
      & all_CoordinateSystem(CoordinateSystemIdx)%CoordinateSystem,Err)

    CALL cmfe_Region_CreateFinish(all_Region(RegionIdx)%Region,Err)

  END DO  ! RegionIdx


  !!!!!!!!!!!!!!!!!!!!!!!!!  Determine if pressure component exist or not in a  Dependent Field      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ALLOCATE(PressureComponentExist(NumberOfDependentField))
  DO DependentFieldIdx = 1, NumberOfDependentField

    IF  ((STR2INT(all_DependentField(DependentFieldIdx)%DependentFieldNumberOfComponents(1)) .GT. &
      MATCH_COORDINATE_SYSTEM(all_CoordinateSystem(DependentFieldIdx)%CoordinateSystemDimension(1)))) THEN  ! if dimensionality of the domain is greater than the number of dependent field components

      PressureComponentExist(DependentFieldIdx) = 1

    ELSE

      PressureComponentExist(DependentFieldIdx) = 0

    END IF

  END DO  !! DependentFieldIdx = 1, NumberOfDependentField
  !!!!!!!!!!!!!!!!!!!!!!!!!  ASSIGNING BASIS TO RESPECTIVE REGIONS       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ALLOCATE(BasisIdxArray(NumberOfGeneratedMesh,2))

  DO GeneratedMeshIdx = 1, NumberOfGeneratedMesh

    CALL MATCH_IDs(all_GeneratedMesh(GeneratedMeshIdx)%GeneratedMeshID(1), &
      & all_DependentField(:)%DependentFieldID(1),DependentFieldIdx, "all_GeneratedMesh(:)", "all_DependentField(:)" )

    CALL MATCH_IDs(all_GeneratedMesh(GeneratedMeshIdx)%GeneratedMeshID(1), &
      & all_Basis(:)%BasisID(1),BasisIdx, "all_GeneratedMesh(:)", "all_Basis(:)" )


    IF (PressureComponentExist(DependentFieldIdx) == 1) THEN

      DO SecondBasisIdx =1, NUmberOfBasis

        IF ((all_Basis(BasisIdx)%BasisID(1) == all_Basis(SecondBasisIdx)%BasisID(1)) .AND. SecondBasisIdx .NE. BasisIdx) THEN

         IF (MATCH_BASIS(all_Basis(BasisIdx)%BasisInterpolationType(1)) .GT.  &
           & MATCH_BASIS(all_Basis(SecondBasisIdx)%BasisInterpolationType(1))) THEN

           BasisIdxArray(GeneratedMeshIdx,1) = BasisIdx

           BasisIdxArray(GeneratedMeshIdx,2) = SecondBasisIdx
        ELSE

           BasisIdxArray(GeneratedMeshIdx,1) = SecondBasisIdx

           BasisIdxArray(GeneratedMeshIdx,2) = BasisIdx

        END IF
      END IF
    END DO

    ELSE

      BasisIdxArray(GeneratedMeshIdx,1) = BasisIdx

    END IF



    DO BasisIdx = 1, PressureComponentExist(DependentFieldIdx) + 1

      j = BasisIdxArray(GeneratedMeshIdx,BasisIdx) !! the real basis index
      
      InterpolationType= MATCH_BASIS(all_Basis(j)%BasisInterpolationType(1))

      NumberGlobalXElements=STR2INT(all_GeneratedMesh(GeneratedMeshIdx)%GeneratedMeshNumberOfElements(1))

      NumberGlobalYElements=STR2INT(all_GeneratedMesh(GeneratedMeshIdx)%GeneratedMeshNumberOfElements(2))

      NumberGlobalZElements=STR2INT(all_GeneratedMesh(GeneratedMeshIdx)%GeneratedMeshNumberOfElements(3))

      !Define geometric basis

      CALL cmfe_Basis_Initialise(all_Basis(j)%Basis,Err)

      CALL cmfe_Basis_CreateStart(BasisUserNumber(j),all_Basis(j)%Basis,Err)

      SELECT CASE(InterpolationType)

        CASE(CMFE_BASIS_LINEAR_LAGRANGE_INTERPOLATION,CMFE_BASIS_QUADRATIC_LAGRANGE_INTERPOLATION, &
          & CMFE_BASIS_CUBIC_LAGRANGE_INTERPOLATION,CMFE_BASIS_CUBIC_HERMITE_INTERPOLATION )

          CALL cmfe_Basis_TypeSet(all_Basis(j)%Basis,CMFE_BASIS_LAGRANGE_HERMITE_TP_TYPE,Err)

        CASE(CMFE_BASIS_QUADRATIC1_HERMITE_INTERPOLATION ,CMFE_BASIS_QUADRATIC2_HERMITE_INTERPOLATION, &
          & CMFE_BASIS_LINEAR_SIMPLEX_INTERPOLATION )

          CALL cmfe_Basis_TypeSet(all_Basis(j)%Basis,CMFE_BASIS_SIMPLEX_TYPE,Err)

      END SELECT

      ! 1D case
      IF(NumberGlobalYElements==0 .AND. NumberGlobalZElements==0) THEN

        CALL cmfe_Basis_NumberOfXiSet(all_Basis(j)%Basis,1,Err)

        CALL cmfe_Basis_InterpolationXiSet(all_Basis(j)%Basis,[InterpolationType],Err)

        CALL cmfe_Basis_QuadratureNumberOfGaussXiSet(all_Basis(j)%Basis, &
          & [STR2INT(all_Basis(j)%BasisNumberOfGaussPoints(1))],Err)

      END IF
       ! 2D case
      IF(NumberGlobalZElements==0) THEN
        CALL cmfe_Basis_NumberOfXiSet(all_Basis(j)%Basis,2,Err)

        CALL cmfe_Basis_InterpolationXiSet(all_Basis(j)%Basis,[InterpolationType,InterpolationType],Err)

        IF (STR2INT(all_Basis(j)%BasisNumberOfGaussPoints(1))>0  &
          .AND. STR2INT(all_Basis(j)%BasisNumberOfGaussPoints(2))>0 ) THEN

          CALL cmfe_Basis_QuadratureNumberOfGaussXiSet(all_Basis(j)%Basis, &
            & [STR2INT(all_Basis(j)%BasisNumberOfGaussPoints(1)), &
              & STR2INT(all_Basis(j)%BasisNumberOfGaussPoints(2))],Err)

        END IF

      ELSE
        CALL cmfe_Basis_NumberOfXiSet(all_Basis(j)%Basis,3,Err)

        CALL cmfe_Basis_InterpolationXiSet(all_Basis(j)%Basis,[InterpolationType,InterpolationType,InterpolationType],Err)

        IF(  STR2INT(all_Basis(j)%BasisNumberOfGaussPoints(1))>0  &
          & .AND.  STR2INT(all_Basis(j)%BasisNumberOfGaussPoints(2))>0  &
            & .AND. STR2INT(all_Basis(j)%BasisNumberOfGaussPoints(3))>0 ) THEN

          CALL cmfe_Basis_QuadratureNumberOfGaussXiSet(all_Basis(j)%Basis, &
            & [STR2INT(all_Basis(j)%BasisNumberOfGaussPoints(1)), &
              &  STR2INT(all_Basis(j)%BasisNumberOfGaussPoints(2)), &
                & STR2INT(all_Basis(j)%BasisNumberOfGaussPoints(3))],Err)

        END IF
      END IF

      CALL cmfe_Basis_CreateFinish(all_Basis(j)%Basis,Err)


    END DO  ! BasisIdx

  END DO !GeneratedMeshIdx
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  DEFINING MESH AND ASSIGNING RESPECTIVE MESHES TO THIER  REGIONS       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  DO GeneratedMeshIdx = 1,NumberOfGeneratedMesh

    !! number of state variables ( for instance 3 ( displacements) +1(Pressure) for 3D cantiliverbeam study).
    NumberOfComponents                   = STR2INT(all_DependentField(GeneratedMeshIdx)%DependentFieldNumberOfComponents(1))

    FirstDimension                       = STR2INT(all_GeneratedMesh(GeneratedMeshIdx)%GeneratedMeshGeometricExtents(1))

    SecondDimension                      = STR2INT(all_GeneratedMesh(GeneratedMeshIdx)%GeneratedMeshGeometricExtents(2))

    ThirdDImension                       = STR2INT(all_GeneratedMesh(GeneratedMeshIdx)%GeneratedMeshGeometricExtents(3))


    !Start the creation of a generated mesh in the region.

    CALL cmfe_GeneratedMesh_Initialise(all_GeneratedMesh(GeneratedMeshIdx)%GeneratedMesh,Err)

    ! MATCH IDS for generated mesh and region objects
    CALL MATCH_IDs(all_GeneratedMesh(GeneratedMeshIdx)%GeneratedMeshID(1), &
       all_Region(:)%RegionID(1),RegionIdx,"all_GeneratedMesh(:)", " all_Region(:)" )

    CALL cmfe_GeneratedMesh_CreateStart(GeneratedMeshUserNumber(GeneratedMeshIdx),all_Region(RegionIdx)%Region, &
     & all_GeneratedMesh(GeneratedMeshIdx)%GeneratedMesh,Err)

    !Set up a regular x*y*z mesh

    CALL cmfe_GeneratedMesh_TypeSet(all_GeneratedMesh(GeneratedMeshIdx)%GeneratedMesh, &
      & MATCH_GENERATED_MESH(all_GeneratedMesh(GeneratedMeshIdx)%GeneratedMeshType(1)),Err)

    ! MATCH objects with similar IDs from the lists: "all_GeneratedMesh(:)" and  "all_DependentField(:)"
    CALL MATCH_IDs(all_GeneratedMesh(GeneratedMeshIdx)%GeneratedMeshID(1), &
       all_DependentField(:)%DependentFieldId(1),DependentFieldIdx,"all_GeneratedMesh(:)", "all_DependentField(:)" )

    IF(PressureComponentExist(DependentFieldIdx)==1) THEN

      !Set the default basis
      CALL cmfe_GeneratedMesh_BasisSet(all_GeneratedMesh(GeneratedMeshIdx)%GeneratedMesh, &
        &  [all_Basis( BasisIdxArray(GeneratedMeshIdx,1))%Basis, &
          & all_Basis( BasisIdxArray(GeneratedMeshIdx,2))%Basis],Err)

    ELSE

      CALL cmfe_GeneratedMesh_BasisSet(all_GeneratedMesh(GeneratedMeshIdx)%GeneratedMesh, &
        & [all_Basis( BasisIdxArray(GeneratedMeshIdx,1))%Basis],Err)

    END IF

    ! 1D case
    IF(NumberGlobalYElements==0 .AND. NumberGlobalZElements==0) THEN

      CALL cmfe_GeneratedMesh_ExtentSet(all_GeneratedMesh(GeneratedMeshIdx)%GeneratedMesh,[FirstDImension],Err)

      CALL cmfe_GeneratedMesh_NumberOfElementsSet(all_GeneratedMesh(GeneratedMeshIdx)%GeneratedMesh,[NumberGlobalXElements],Err)

    ! 2D case
    ELSE IF (NumberGlobalZElements==0) THEN

      CALL cmfe_GeneratedMesh_ExtentSet(all_GeneratedMesh(GeneratedMeshIdx)%GeneratedMesh,[FirstDImension,SecondDimension],Err)

      CALL cmfe_GeneratedMesh_OriginSet(all_GeneratedMesh(GeneratedMeshIdx)%GeneratedMesh, &
        & [STR2REAL(all_GeneratedMesh(GeneratedMeshIdx)%GeneratedMeshOriginSet(1)), &
          & STR2REAL(all_GeneratedMesh(GeneratedMeshIdx)%GeneratedMeshOriginSet(2))],Err)

      CALL cmfe_GeneratedMesh_NumberOfElementsSet(all_GeneratedMesh(GeneratedMeshIdx)%GeneratedMesh, &
        & [NumberGlobalXElements,NumberGlobalYElements],Err)

    ELSE

      CALL cmfe_GeneratedMesh_ExtentSet(all_GeneratedMesh(GeneratedMeshIdx)%GeneratedMesh, &
        & [FirstDImension,SecondDimension,ThirdDImension],Err)

      CALL cmfe_GeneratedMesh_OriginSet(all_GeneratedMesh(GeneratedMeshIdx)%GeneratedMesh, &
        & [STR2REAL(all_GeneratedMesh(GeneratedMeshIdx)%GeneratedMeshOriginSet(1)), &
          & STR2REAL(all_GeneratedMesh(GeneratedMeshIdx)%GeneratedMeshOriginSet(2)), &
            & STR2REAL(all_GeneratedMesh(GeneratedMeshIdx)%GeneratedMeshOriginSet(3))],Err)

      CALL cmfe_GeneratedMesh_NumberOfElementsSet(all_GeneratedMesh(GeneratedMeshIdx)%GeneratedMesh, &
       & [NumberGlobalXElements,NumberGlobalYElements, NumberGlobalZElements],Err)

    END IF

    ! Match IDs of the GeneratedMesh and Mesh object
    CALL MATCH_IDs(all_GeneratedMesh(GeneratedMeshIdx)%GeneratedMeshID(1), &
       all_Mesh(:)%MeshId(1),MeshIdx,"all_GeneratedMesh(:)", "all_Mesh(:)" )

    CALL cmfe_Mesh_Initialise(all_Mesh(MeshIdx)%Mesh,Err)

    !Finish the creation of a generated mesh in the region
    CALL cmfe_GeneratedMesh_CreateFinish(all_GeneratedMesh(GeneratedMeshIdx)%GeneratedMesh,MeshUserNumber(MeshIdx), &
      & all_Mesh(MeshIdx)%Mesh,Err)

  END DO !GeneratedMeshIdx

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  DECOMPOSING BLOCK           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  DO DecompositionIdx = 1,NumberOfDecomposition
    !Create a decomposition

    CALL cmfe_Decomposition_Initialise(all_Decomposition(DecompositionIdx)%Decomposition,Err)

    CALL cmfe_Decomposition_CreateStart(DecompositionUserNumber(DecompositionIdx), &
      & all_Mesh(DecompositionIdx)%Mesh,all_Decomposition(DecompositionIdx)%Decomposition,Err)

    CALL cmfe_Decomposition_TypeSet(all_Decomposition(DecompositionIdx)%Decomposition,CMFE_DECOMPOSITION_CALCULATED_TYPE,Err)

    CALL cmfe_Decomposition_NumberOfDomainsSet(all_Decomposition(DecompositionIdx)%Decomposition,NumberOfDomains,Err)

    IF (all_Decomposition(NumberOfDecomposition)%CalculateElementFaces(1) == "TRUE") THEN

      CALL cmfe_Decomposition_CalculateFacesSet(all_Decomposition(DecompositionIdx)%Decomposition,.TRUE.,Err)

    END IF

    CALL cmfe_Decomposition_CreateFinish(all_Decomposition(DecompositionIdx)%Decomposition,Err)

  END DO !DecompositionIdx
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   GEOMETRIC  BLOCK           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  DO GeometricFieldIdx  = 1, NumberOfGeometricField
    
    ! Match IDs of the GeometricField and Decomposition object

    CALL MATCH_IDs(all_GeometricField(GeometricFieldIdx)%GeometricFieldID(1), &
       all_Decomposition(:)%DecompositionID(1),DecompositionIdx,"all_GeometricField", "all_Decomposition(:)" )

    ! Match IDs of the GeometricField and Region object

    CALL MATCH_IDs(all_GeometricField(GeometricFieldIdx)%GeometricFieldID(1), &
      & all_Region(:)%RegionID(1),RegionIdx,"all_GeometricField", "all_Region(:)" )

    ! Match IDs of the GeometricField and GeneratedMesh object

    CALL MATCH_IDs(all_GeometricField(GeometricFieldIdx)%GeometricFieldID(1), &
       all_GeneratedMesh(:)%GeneratedMeshID(1),GeneratedMeshIdx,"all_GeometricField(:)", "all_GeneratedMesh(:)" )

    !Create a field to put the geometry (default is geometry)

    CALL cmfe_Field_Initialise(all_GeometricField(GeometricFieldIdx)%GeometricField,Err)

    CALL cmfe_Field_CreateStart(FieldGeometryUserNumber(GeometricFieldIdx),all_Region(RegionIdx)%Region, &
      & all_GeometricField(GeometricFieldIdx)%GeometricField,Err)

    CALL cmfe_Field_MeshDecompositionSet(all_GeometricField(GeometricFieldIdx)%GeometricField, &
      & all_Decomposition(DecompositionIdx)%Decomposition,Err)

    CALL cmfe_Field_VariableLabelSet(all_GeometricField(GeometricFieldIdx)%GeometricField, &
      & CMFE_FIELD_U_VARIABLE_TYPE,all_GeometricField(GeometricFieldIdx)%GeometricFieldLabel(1),Err)

    CALL cmfe_Field_ScalingTypeSet(all_GeometricField(GeometricFieldIdx)%GeometricField, &
      & CMFE_FIELD_ARITHMETIC_MEAN_SCALING,Err)

    CALL cmfe_Field_CreateFinish(all_GeometricField(GeometricFieldIdx)%GeometricField,Err)



    !Update the geometric field parameters
    CALL cmfe_GeneratedMesh_GeometricParametersCalculate(all_GeneratedMesh(GeneratedMeshIdx)%GeneratedMesh, &
      & all_GeometricField(GeometricFieldIdx)%GeometricField,Err)
  END DO ! GeometricFieldIdx  !!!!! END LOOP FOR GEOMETRIC FIELD

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   FIBER FIELD  BLOCK           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  DO FiberFieldIdx = 1, NumberOfFiberField


    ! Match fiber field objects with similar ID objects  of other data type
    CALL MATCH_IDs(all_FibreField(FiberFieldIdx)%FibreFieldID(1), &
       all_Region(:)%RegionID(1),RegionIdx,"all_FiberField(:)", "all_Region(:)" )

    CALL MATCH_IDs(all_FibreField(FiberFieldIdx)%FibreFieldID(1), &
       all_Decomposition(:)%DecompositionID(1),DecompositionIdx,"all_FiberField(:)", "all_Decomposition(:)" )

    CALL MATCH_IDs(all_FibreField(FiberFieldIdx)%FibreFieldID(1), &
       all_GeometricField(:)%GeometricFieldID(1) ,GeometricFieldIdx, "all_FiberField(:)", "all_GeometricField(:)" )

    CALL cmfe_Field_Initialise(all_FibreField(FiberFieldIdx)%FibreField,Err)

    CALL cmfe_Field_CreateStart(FieldFibreUserNumber(FiberFieldIdx),all_Region(RegionIdx)%Region, &
      & all_FibreField(FiberFieldIdx)%FibreField,Err)

    CALL cmfe_Field_TypeSet(all_FibreField(FiberFieldIdx)%FibreField,CMFE_FIELD_FIBRE_TYPE,Err)

    CALL cmfe_Field_MeshDecompositionSet(all_FibreField(FiberFieldIdx)%FibreField,  &
      all_Decomposition(DecompositionIdx)%Decomposition,Err)

    CALL cmfe_Field_GeometricFieldSet(all_FibreField(FiberFieldIdx)%FibreField, &
      all_GeometricField(GeometricFieldIdx)%GeometricField,Err)

    CALL cmfe_Field_VariableLabelSet(all_FibreField(FiberFieldIdx)%FibreField,CMFE_FIELD_U_VARIABLE_TYPE, &
      & all_FibreField(FiberFieldIdx)%FiberFieldLabel(1),Err)

    CALL cmfe_Field_ScalingTypeSet(all_FibreField(FiberFieldIdx)%FibreField,CMFE_FIELD_ARITHMETIC_MEAN_SCALING,Err)
    CALL cmfe_Field_CreateFinish(all_FibreField(FiberFieldIdx)%FibreField,Err)

    !Set the fibre field

    DO ComponentIdx=1,NumberOfComponents - PressureComponentExist(RegionIdx) 
      CALL cmfe_Field_ComponentValuesInitialise(all_FibreField(FiberFieldIdx)%FibreField, &
        & CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1, &
          & STR2REAL(all_FibreField(FiberFieldIdx)%FibreFieldParameters(ComponentIdx)),Err)
    END DO ! ComponentIdx

  END DO   !!!!! END LOOP FOR Fiber FIELD

  DO EquationSetIdx  = 1,NumberOfEquationsSet            !!!! loop for eqaution set

    CALL MATCH_IDs(all_EquationsSet(EquationSetIdx)%EquationSetID(1), &
       all_Region(:)%RegionID(1),RegionIdx,"all_EquationsSet(:)", "all_Region(:)" )

    CALL cmfe_EquationsSet_Initialise(all_EquationsSet(EquationSetIdx)%EquationsSet,Err)

    CALL cmfe_Field_Initialise(all_EquationsSetField(EquationSetIdx)%EquationsSetField,Err)

    IF (NumberOfFiberField .GE. 1) THEN

      CALL MATCH_IDs(all_EquationsSet(EquationSetIdx)%EquationSetID(1), &
        all_FibreField(:)%FibreFieldID(1),FiberFieldIdx,"all_EquationsSet(:)", "all_Fiber(:)" )

      CALL cmfe_EquationsSet_CreateStart(EquationSetUserNumber(EquationSetIdx),all_Region(RegionIdx)%Region, &
        & all_FibreField(FiberFieldIdx)%FibreField, &
          & [MATCH_EQUATION_SET(all_EquationsSet(EquationSetIdx)%EquationSetClass(1)), &
            & MATCH_EQUATION_SET(all_EquationsSet(EquationSetIdx)%EquationSetType(1)), &
              & MATCH_EQUATION_SET(all_EquationsSet(EquationSetIdx)%EquationSetSubType(1))], &
                & EquationsSetFieldUserNumber(EquationSetIdx),all_EquationsSetField(EquationSetIdx)%EquationsSetField, &
                  & all_EquationsSet(EquationSetIdx)%EquationsSet,Err)


     ELSE

      CALL MATCH_IDs(all_EquationsSet(EquationSetIdx)%EquationSetID(1), &
        all_GeometricField(:)%GeometricFieldID(1),GeometricFieldIdx,"all_EquationsSet(:)", "all_GeometricField(:)" )

      CALL cmfe_EquationsSet_CreateStart(EquationSetUserNumber(EquationSetIdx),all_Region(RegionIdx)%Region, &
        & all_GeometricField(GeometricFieldIdx)%GeometricField, &
          & [MATCH_EQUATION_SET(all_EquationsSet(EquationSetIdx)%EquationSetClass(1)), &
            & MATCH_EQUATION_SET(all_EquationsSet(EquationSetIdx)%EquationSetType(1)), &
              & MATCH_EQUATION_SET(all_EquationsSet(EquationSetIdx)%EquationSetSubType(1))], &
                & EquationsSetFieldUserNumber(EquationSetIdx),all_EquationsSetField(EquationSetIdx)%EquationsSetField, &
                  & all_EquationsSet(EquationSetIdx)%EquationsSet,Err)


    END IF

    CALL cmfe_EquationsSet_CreateFinish(all_EquationsSet(EquationSetIdx)%EquationsSet,Err)

  END DO ! EquationSetIdx


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   Dependent  FIELD  BLOCK           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !Create the dependent field
  DO DependentFieldIdx = 1, NumberOfDependentField
    !Match ID of DependentField and EquationSet objects

    CALL MATCH_IDs(all_DependentField(DependentFieldIdx)%DependentFieldID(1), &
       all_EquationsSet(:)%EquationSetID(1),EquationSetIdx,"all_DependentField(:)", "all_EquationSet(:)" )


    !Create the Dependent field
    CALL cmfe_Field_Initialise(all_DependentField(DependentFieldIdx)%DependentField,Err)

    CALL cmfe_EquationsSet_DependentCreateStart(all_EquationsSet(EquationSetIdx)%EquationsSet, &
      & FieldDependentUserNumber(DependentFieldIdx), all_DependentField(DependentFieldIdx)%DependentField,Err)

    CALL cmfe_Field_VariableLabelSet(all_DependentField(DependentFieldIdx)%DependentField, &
      & CMFE_FIELD_U_VARIABLE_TYPE,all_DependentField(DependentFieldIdx)%DependentFieldLabel(1),Err)

    DO ComponentIdx=1,NumberOfComponents - PressureComponentExist(DependentFieldIdx)   !! "-PressureComponent" to exclude pressure from the list of components 

      CALL cmfe_Field_ComponentMeshComponentSet(all_DependentField(DependentFieldIdx)%DependentField, &
        & CMFE_FIELD_U_VARIABLE_TYPE,ComponentIdx,1,Err)

      CALL cmfe_Field_ComponentMeshComponentSet(all_DependentField(DependentFieldIdx)%DependentField, &
        & CMFE_FIELD_DELUDELN_VARIABLE_TYPE,ComponentIdx,1,Err)

    END DO

    IF(NumberofBasis == 2) THEN           !! HARD CODED
      IF (NumberOfComponents .GT. 3) THEN !! HARD CODED DUE TO AN ISSUE, WILL BE REMOVED AFTER THE NEXT MEETING.
        CALL cmfe_Field_ComponentInterpolationSet(all_DependentField(DependentFieldIdx)%DependentField, &
          & CMFE_FIELD_U_VARIABLE_TYPE,NumberOfComponents,CMFE_FIELD_NODE_BASED_INTERPOLATION,Err)
        CALL cmfe_Field_ComponentInterpolationSet(all_DependentField(DependentFieldIdx)%DependentField, &
          & CMFE_FIELD_DELUDELN_VARIABLE_TYPE,NumberOfComponents, CMFE_FIELD_NODE_BASED_INTERPOLATION,Err)
      END IF

      CALL cmfe_Field_ComponentMeshComponentSet( &
        & all_DependentField(DependentFieldIdx)%DependentField,CMFE_FIELD_U_VARIABLE_TYPE, &
          & NumberOfComponents,2,Err)
      CALL cmfe_Field_ComponentMeshComponentSet( &
        & all_DependentField(DependentFieldIdx)%DependentField,CMFE_FIELD_DELUDELN_VARIABLE_TYPE, &
          & NumberOfComponents,2,Err)

    END IF

    CALL cmfe_EquationsSet_DependentCreateFinish(all_EquationsSet(DependentFieldIdx)%EquationsSet,Err)

    IF  (TRIM(all_DependentField(DependentFieldIdx)%DependentFieldInitialValueOfStateVector(1)) == "UNDEFORMED") THEN

      DO ComponentIdx=1,NumberOfComponents - PressureComponentExist(DependentFieldIdx)
        CALL cmfe_Field_ParametersToFieldParametersComponentCopy(all_GeometricField(DependentFieldIdx)%GeometricField, &
          & CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
            & ComponentIdx,all_DependentField(DependentFieldIdx)%DependentField, &
              & CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,ComponentIdx,Err)
      END DO


    ELSE
      DO ComponentIdx=1,NumberOfComponents - PressureComponentExist(DependentFieldIdx)

        CALL cmfe_Field_ComponentValuesInitialise(all_DependentField(DependentFieldIdx)%DependentField, &
          &  CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,ComponentIdx, &
            & STR2REAL(all_DependentField(DependentFieldIdx)%DependentFieldInitialValueOfStateVector(ComponentIdx)),Err)

      END DO  ! ComponentIdx
    END IF  ! TRIM(DependentField_arguments(4,DependentFieldIdx)) == "UNDEFORMED")

    IF(PressureComponentExist(DependentFieldIdx) == 1) THEN

      CALL cmfe_Field_ComponentValuesInitialise(all_DependentField(DependentFieldIdx)%DependentField, &
        &  CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,NUmberOfComponents, &
          & STR2REAL(all_DependentField(DependentFieldIdx)%DependentFieldInitialValueOfStateScalar(1)),Err)

    END IF

    CALL cmfe_Field_ParameterSetUpdateStart(all_DependentField(DependentFieldIdx)%DependentField,CMFE_FIELD_U_VARIABLE_TYPE &
      & ,CMFE_FIELD_VALUES_SET_TYPE,Err)
    CALL cmfe_Field_ParameterSetUpdateFinish(all_DependentField(DependentFieldIdx)%DependentField,CMFE_FIELD_U_VARIABLE_TYPE, &
      & CMFE_FIELD_VALUES_SET_TYPE,Err)

  END DO ! DependentFieldIdxx

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   MATERIAL FIELD  BLOCK           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  DO MaterialFieldIdx  = 1,NumberOfMaterialField

    CALL MATCH_IDs(all_MaterialField(MaterialFieldIdx)%MaterialFieldID(1), &
       all_EquationsSet(:)%EquationSetID(1),EquationSetIdx,"all_MaterialField(:)", "all_EquationSet(:)" )

    CALL cmfe_Field_Initialise(all_MaterialField(EquationSetIdx)%MaterialField,Err)

    CALL cmfe_EquationsSet_MaterialsCreateStart(all_EquationsSet(EquationSetIdx)%EquationsSet, &
        & FieldMaterialUserNumber(EquationSetIdx),all_MaterialField(MaterialFieldIdx)%MaterialField,Err)

    CALL cmfe_Field_VariableLabelSet(all_MaterialField(MaterialFieldIdx)%MaterialField, &
      & CMFE_FIELD_U_VARIABLE_TYPE,all_MaterialField(MaterialFieldIdx)%MaterialFieldLabel(1),Err)

    CALL cmfe_EquationsSet_MaterialsCreateFinish(all_EquationsSet(EquationSetIdx)%EquationsSet,Err)

    DO  MaterialParametersIdx = 1,MATCH_MATERIAL_PARAMETERS(all_EquationsSet(EquationSetIdx)%EquationSetSubType(1))

      CALL  cmfe_Field_ComponentValuesInitialise(all_MaterialField(MaterialFieldIdx)%MaterialField, &
        & CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, MaterialParametersIdx , &
          & STR2REAL(all_MaterialField(MaterialFieldIdx)%MaterialFieldParameters(MaterialParametersIdx)),Err)
    END DO

  END DO  ! MaterialFieldIdx

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   EQUATION SET BLOCK           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  DO EquationSetIdx  = 1,NumberOfEquationsSet            !!!! loop for eqaution set

    !Create the equations set equations
    CALL cmfe_Equations_Initialise(all_Equations%Equations(EquationSetIdx),Err)

    CALL cmfe_EquationsSet_EquationsCreateStart(all_EquationsSet(EquationSetIdx)%EquationsSet, &
      & all_Equations%Equations(EquationSetIdx),Err)

    CALL cmfe_Equations_SparsityTypeSet(all_Equations%Equations(EquationSetIdx),CMFE_EQUATIONS_SPARSE_MATRICES,Err)

    CALL cmfe_Equations_OutputTypeSet(all_Equations%Equations(EquationSetIdx), &
      & MATCH_OUTPUT(all_EquationsSet(EquationSetIdx)%EquationSetOutputTypeSet(1)),Err)

    CALL cmfe_EquationsSet_EquationsCreateFinish(all_EquationsSet(EquationSetIdx)%EquationsSet,Err)


    IF (NumberOfSourceField .GT. 0) THEN
      print *, all_SourceField(1)%SourceFieldID(1)
      CALL MATCH_IDs(all_EquationsSet(EquationSetIdx)%EquationSetID(1), &
        & all_SourceField(:)%SourceFieldID(1),SourceFieldIdx,"all_EquationsSet(:)", "all_Source(:)" )
      IF (all_SourceField(SourceFieldIdx)%SourceFieldType(1) == "GRAVITY" ) THEN   !! look input file IF the given field is gravitational
      !Create the source field with the gravi    print *, "hello" , CoordinateSystem(EquationSetIdx)%CoordinateSystemDimension(1)ty vector


        CALL cmfe_Field_Initialise(all_SourceField(SourceFieldIdx)%SourceField,Err)

        CALL cmfe_EquationsSet_SourceCreateStart(all_EquationsSet(EquationSetIdx)%EquationsSet, &
            & FieldUserNumber(NumberOfSourceField) ,all_SourceField(SourceFieldIdx)%SourceField,Err)

        CALL cmfe_Field_ScalingTypeSet(all_SourceField(SourceFieldIdx)%SourceField,CMFE_FIELD_ARITHMETIC_MEAN_SCALING,Err)

        CALL cmfe_EquationsSet_SourceCreateFinish(all_EquationsSet(EquationSetIdx)%EquationsSet,Err)

        DO ComponentIdx=1,MATCH_COORDINATE_SYSTEM(all_CoordinateSystem(CoordinateSystemIdx)%CoordinateSystemDimension(1))  !!  MATCH_COORDINATE_SYSTEM(CoordinateSystem(CoordinateSystemIdx)%
                                                                                                                           !!  ...CoordinateSystemDimension(1)) will give 2 for 2D structure and 3 for 3D
                                                                                                                           !!  ....structure

          CALL cmfe_Field_ComponentValuesInitialise(all_SourceField(SourceFieldIdx)%SourceField,  &
            & CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,ComponentIdx, &
              & STR2REAL(all_SourceField(SourceFieldIdx)%SourceFieldComponents(ComponentIdx)),Err)

        END DO ! ComponentIdx

      END IF
    ENDIF !(NumberOfSourceField .GT. 0)
  END DO ! EquationsSetIdx
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   PROBLEM BLOCK           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  DO ProblemIdx = 1, NumberOfProblem

    !Define the problem

    CALL cmfe_Problem_Initialise(all_Problem(ProblemIdx)%Problem,Err)

    CALL cmfe_Problem_CreateStart(ProblemUserNumber(ProblemIdx),[MATCH_PROBLEM(all_Problem(ProblemIdx)%ProblemClass(1)), &
      & MATCH_PROBLEM(all_Problem(ProblemIdx)%ProblemType(1)),MATCH_PROBLEM(all_Problem(ProblemIdx)%ProblemSubType(1))], &
        & all_Problem(ProblemIdx)%Problem,Err)

    CALL cmfe_Problem_CreateFinish(all_Problem(ProblemIdx)%Problem,Err)

  END DO ! ProblemIdx
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   CONTROL LOOP BLOCK           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  DO ControlLoopIdx = 1,NumberOfControlLoop
    !Create the problem control loop


    CALL MATCH_IDs(all_ControlLoop(ControlLoopIdx)%ControlLoopID(1), &
      & all_Problem(:)%ProblemID(1),ProblemIdx,"all_ControlLoop()", "all_Problem()" )

    CALL cmfe_Problem_ControlLoopCreateStart(all_Problem(ProblemIdx)%Problem,Err)

    CALL cmfe_ControlLoop_Initialise(all_ControlLoop(ControlLoopIdx)%ControlLoop,Err)

    CALL cmfe_Problem_ControlLoopGet(all_Problem(ProblemIdx)%Problem,CMFE_CONTROL_LOOP_NODE, &
      & all_ControlLoop(ControlLoopIdx)%ControlLoop,Err)

    CALL cmfe_ControlLoop_TypeSet(all_ControlLoop(ControlLoopIdx)%ControlLoop, &
      & MATCH_CONTROL_LOOP(all_ControlLoop(ControlLoopIdx)%ControlLoopType(1)),Err)

    SELECT CASE(MATCH_CONTROL_LOOP(all_ControlLoop(ControlLoopIdx)%ControlLoopType(1)))

      CASE (CMFE_PROBLEM_CONTROL_LOAD_INCREMENT_LOOP_TYPE)

        CALL cmfe_ControlLoop_MaximumIterationsSet(all_ControlLoop(ControlLoopIdx)%ControlLoop, &
          & STR2INT(all_ControlLoop(ControlLoopIdx)%ControlLoopLoadIncrement(1)),Err)

        CALL cmfe_ControlLoop_LoaDOutputSet(all_ControlLoop(ControlLoopIdx)%ControlLoop,1,Err)

      CASE (CMFE_PROBLEM_CONTROL_TIME_LOOP_TYPE)

        CALL cmfe_ControlLoop_TimesSet(all_ControlLoop(ControlLoopIdx)%ControlLoop, &
          & STR2REAL(all_ControlLoop(ControlLoopIdx)%ControlLoopTimeIncrement(1)), &
            & STR2REAL(all_ControlLoop(ControlLoopIdx)%ControlLoopTimeIncrement(2)), &
              &  STR2REAL(all_ControlLoop(ControlLoopIdx)%ControlLoopTimeIncrement(3)),Err)

        CALL cmfe_ControlLoop_TimeOutputSet(all_ControlLoop(ControlLoopIdx)%ControlLoop,1,Err)

    END SELECT

    CALL cmfe_Problem_ControlLoopCreateFinish(all_Problem(ProblemIdx)%Problem,Err)

  END DO ! ControlLoopIdx



  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   SOLVER BLOCK           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  DO solverIdx = 1,NumberOfSolver

    !Create the problem solvers
    CALL cmfe_Solver_Initialise(all_Solver(solverIdx)%Solver,Err)
    CALL cmfe_Solver_Initialise(all_LinearSolver(solverIdx)%LinearSolver,Err)
    CALL cmfe_Problem_SolversCreateStart(all_Problem(solverIdx)%Problem,Err)
    CALL cmfe_Problem_SolverGet(all_Problem(solverIdx)%Problem,CMFE_CONTROL_LOOP_NODE,1,all_Solver(solverIdx)%Solver,Err)
    CALL cmfe_Solver_OutputTypeSet(all_Solver(solverIdx)%Solver,MATCH_OUTPUT(all_Solver(solverIdx)%OutputTypeSet(1)),Err)
    IF ((TRIM(all_Solver(solverIdx)%NonlinearSolverType(1)) == "STATIC_NONLINEAR") .OR. &
      & (TRIM(all_Solver(solverIdx)%NonlinearSolverType(1)) == "DYNAMIC_NONLINEAR") ) THEN

      CALL cmfe_Solver_Initialise(all_NonLinearSolver(solverIdx)%NonLinearSolver,Err)

      IF (TRIM(all_Solver(solverIdx)%NonlinearSolverType(1)) == "STATIC_NONLINEAR") THEN
        all_NonLinearSolver(solverIdx)%NonLinearSolver = all_Solver(solverIdx)%Solver

      ELSE IF (TRIM(all_Solver(solverIdx)%NonlinearSolverType(1)) == "DYNAMIC_NONLINEAR") THEN
        !! using implicit /theeta scheme for /theeta = 1
        CALL cmfe_Solver_DynamicThetaSet(all_Solver(solverIdx)%Solver,1.0_CMISSRP,Err)

        CALL cmfe_Solver_DynamicNonlinearSolverGet(all_Solver(solverIdx)%Solver,all_NonlinearSolver(solverIdx)%NonlinearSolver,Err)

      END IF

      CALL cmfe_Solver_NewtonJacobianCalculationTypeSet(all_NonLinearSolver(solverIdx)%NonLinearSolver, &
        & MATCH_SOLVER(all_Solver(solverIdx)%JacobianType(1)),Err)


      CALL cmfe_Solver_NewtonLinearSolverGet(all_NonlinearSolver(solverIdx)%NonlinearSolver, &
        & all_LinearSolver(solverIdx)%LinearSolver,Err)
      
      IF (STR2REAL(all_Solver(solverIdx)%NewtonRelativeTolerance(1)) .NE. 0) THEN !! default value

        CALL cmfe_Solver_NewtonRelativeToleranceSet(all_NonlinearSolver(solverIdx)%NonlinearSolver, &
          & STR2REAL(all_Solver(solverIdx)%NewtonRelativeTolerance(1)),Err)

      END IF

      IF (STR2REAL(all_Solver(solverIdx)%NewtonAbsoluteTolerance(1)) .NE. 0) THEN !! default value

      CALL cmfe_Solver_NewtonAbsoluteToleranceSet(all_NonlinearSolver(solverIdx)%NonlinearSolver, &
        & STR2REAL(all_Solver(solverIdx)%NewtonAbsoluteTolerance(1)),Err)

      END IF

      CALL cmfe_Solver_NewtonMaximumIterationsSet(all_NonlinearSolver(solverIdx)%NonlinearSolver, &
        & STR2INT(all_Solver(solverIdx)%NewtonMaximumIterations(1)),Err)

      CALL cmfe_Solver_OutputTypeSet(all_NonlinearSolver(solverIdx)%NonlinearSolver, &
        & MATCH_OUTPUT(all_Solver(solverIdx)%OutputTypeSet(1)),Err)

      CALL cmfe_Solver_OutputTypeSet(all_LinearSolver(solverIdx)%LinearSolver, &
        & MATCH_OUTPUT(all_Solver(solverIdx)%OutputTypeSet(1)),Err)

      !! IF THE SPECIFIED LINEAER SOLVER IS DIRECT

      IF (MATCH_SOLVER(all_Solver(solverIdx)%LinearSolverType(1)) == CMFE_SOLVER_LINEAR_DIRECT_SOLVE_TYPE ) THEN

        CALL cmfe_Solver_LinearTypeSet(all_LinearSolver(solverIdx)%LinearSolver, &
          & MATCH_SOLVER(all_Solver(solverIdx)%LinearSolverType(1)),Err)

        CALL cmfe_Solver_LinearDirectTypeSet(all_LinearSolver(solverIdx)%LinearSolver, &
          & MATCH_SOLVER(all_Solver(solverIdx)%DirectSolverType(1)),Err)

        CALL cmfe_Solver_LibraryTypeSet(all_LinearSolver(solverIdx)%LinearSolver, &
          & MATCH_SOLVER(all_Solver(solverIdx)%LibraryType(1)),Err)

      !! IF THE SPECIFIED LINEAER SOLVER IS ITERATIVE
      ELSE IF  (MATCH_SOLVER(all_Solver(solverIdx)%LinearSolverType(1)) == CMFE_SOLVER_LINEAR_ITERATIVE_SOLVE_TYPE ) THEN

        CALL cmfe_Solver_LinearTypeSet(all_LinearSolver(solverIdx)%LinearSolver, &
          & MATCH_SOLVER(all_Solver(solverIdx)%LinearSolverType(1)),Err)

        CALL cmfe_Solver_LinearIterativeTypeSet(all_LinearSolver(solverIdx)%LinearSolver, &
          & MATCH_SOLVER(all_Solver(solverIdx)%IterativeSolverType(1)),Err)

        IF (all_Solver(solverIdx)%IterativeSolverAbsoluteTolerance(1) == "0" ) THEN !! default value 

        CALL cmfe_Solver_LinearIterativeAbsoluteToleranceSet(all_LinearSolver(solverIdx)%LinearSolver, &
          & STR2REAL(all_Solver(solverIdx)%IterativeSolverAbsoluteTolerance(1)),Err)

        END IF

        IF (all_Solver(solverIdx)%IterativeSolverRelativeTolerance(1) == "0" ) THEN !! default value

          CALL cmfe_Solver_LinearIterativeRelativeToleranceSet(all_LinearSolver(solverIdx)%LinearSolver, &
            & STR2REAL(all_Solver(solverIdx)%IterativeSolverRelativeTolerance(1)),Err)

        END IF

        CALL cmfe_Solver_LinearIterativeMaximumIterationsSet(all_LinearSolver(solverIdx)%LinearSolver, &
          & STR2INT(all_Solver(solverIdx)%IterativeSolverMaximumIterations(1)),Err)

        CALL cmfe_Solver_LinearIterativePreconditionerTypeSet(all_LinearSolver(solverIdx)%LinearSolver, &
          & MATCH_SOLVER(all_Solver(solverIdx)%IterativePreconditioner(1)) ,Err)

      END IF


    END IF ! IF ((TRIM(Solver%NonlinearSolver(1,solverIdx)) == "STATIC_NONLINEAR") .OR. &

    CALL cmfe_Problem_SolversCreateFinish(all_Problem(solverIdx)%Problem,Err)

    !Create the problem solver equations
    CALL cmfe_Solver_Initialise(all_Solver(solverIdx)%Solver,Err)
    CALL cmfe_SolverEquations_Initialise(all_SolverEquations(solverIdx)%SolverEquations,Err)
    CALL cmfe_Problem_SolverEquationsCreateStart(all_Problem(solverIdx)%Problem,Err)
    CALL cmfe_Problem_SolverGet(all_Problem(solverIdx)%Problem,CMFE_CONTROL_LOOP_NODE,1,all_Solver(solverIdx)%Solver,Err)
    CALL cmfe_Solver_SolverEquationsGet(all_Solver(solverIdx)%Solver,all_SolverEquations(solverIdx)%SolverEquations,Err)
    CALL cmfe_SolverEquations_EquationsSetAdd(all_SolverEquations(solverIdx)%SolverEquations, &
                                               & all_EquationsSet(solverIdx)%EquationsSet, EquationsSetIndex,Err)
    CALL cmfe_Problem_SolverEquationsCreateFinish(all_Problem(solverIdx)%Problem,Err)

  END DO ! solverIdx 


  DO BOundaryConditionIdx = 1,NumberOfBoundaryCondition

    !Prescribe boundary conditions (absolute nodal parameters)
    CALL cmfe_BoundaryConditions_Initialise(all_BoundaryConditions(BOundaryConditionIdx)%BoundaryConditions,Err)
    CALL cmfe_SolverEquations_BoundaryConditionsCreateStart(all_SolverEquations(BOundaryConditionIdx)%SolverEquations, &
       & all_BoundaryConditions(BOundaryConditionIdx)%BoundaryConditions,Err)

    DO i = 1,NumberOfDirichelet

      CALL cmfe_GeneratedMesh_SurfaceGet(all_GeneratedMesh(1)%GeneratedMesh, &
        & MATCH_BC(BoundaryDirichelet(2+6*(i-1),BOundaryConditionIdx)),SurfaceNodes,LeftNormalXi,Err)
      if (i == 5) THEN  !!! the following lines are hard coded and will be removed later upon discussion with Andreas
        h = 2
      else
        h = 1
      end if

      DO nodeIdx=1,SIZE(SurfaceNodes,1),h
      NodeNumber=SurfaceNodes(nodeIdx)
      !! extract spatial coordinates for NodeNumber
      CALL cmfe_Field_ParameterSetGetNode(all_GeometricField(1)%GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,1, &
        & 1,1,NodeNumber,1,x,Err)
      
    ! 1D case
      IF(NumberGlobalYElements==0 .AND. NumberGlobalZElements==0) THEN
        y = 0
        z = 0
      ELSE
        CALL cmfe_Field_ParameterSetGetNode(all_GeometricField(1)%GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,1, &
          & 1,1,NodeNumber,2,y,Err)
      END IF

    ! 2D case
      IF (NumberGlobalZElements==0) THEN
        z = 0
      ELSE
        CALL cmfe_Field_ParameterSetGetNode(all_GeometricField(1)%GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,1, &
          & 1,1,NodeNumber,3,z,Err)      
      END IF

      FunctionName = BoundaryDirichelet(4+6*(i-1),BOundaryConditionIdx)
      IF (FunctionName(1:9) == "QUADRATIC") THEN

        !! FunctionName is the function ID
        !! FInd FUnction with label FunctionName in list of functions all_Function(:)
        CALL MATCH_IDs(TRIM(FunctionName),all_Function(:)%FUnctionID(1),FunctionIdx)
        CALL QUADRATIC(x,y,z, STR2REAL(all_Function(FunctionIdx)%FunctionConstants(1)), &
          & STR2REAL(all_Function(FunctionIdx)%FunctionConstants(2)), &
            & STR2REAL(all_Function(FunctionIdx)%FunctionConstants(3)), &
              &  STR2REAL(all_Function(FunctionIdx)%FunctionConstants(4)), &
                & STR2REAL(all_Function(FunctionIdx)%FunctionConstants(5)), &
                  & STR2REAL(all_Function(FunctionIdx)%FunctionConstants(6)),&
                    & STR2REAL(all_Function(FunctionIdx)%FunctionConstants(7)), &
                      & STR2REAL(all_Function(FunctionIdx)%FunctionConstants(8)), &
                        & STR2REAL(all_Function(FunctionIdx)%FunctionConstants(9)), &
                          & STR2REAL(all_Function(FunctionIdx)%FunctionConstants(10)), &
                            & BOundaryConditionValue)
      ELSE

        BOundaryConditionValue = STR2REAL(BoundaryDirichelet(4+6*(i-1),BOundaryConditionIdx)) !! for uniform BC
      END IF


      !!!!!!!!

      CALL cmfe_Decomposition_NodeDOmainGet(all_Decomposition(1)%Decomposition,NodeNumber,1,NodeDOmain,Err)
      IF(NodeDOmain==ComputationalNodeNumber) THEN

        DO ComponentIdx=1,NumberOfComponents
          Constraint =  BoundaryDirichelet(3+6*(i-1),BoundaryConditionIdx)

          IF (TRIM(Constraint(ComponentIdx:ComponentIdx)) == "1") THEN
            ! decide whether we add to the nodal value or set it to a given value
            IF (TRIM(BoundaryDirichelet(6+6*(i-1),BOundaryConditionIdx)) == "SET") THEN
              CALL cmfe_BoundaryConditions_SetNode( &
                & all_BoundaryConditions(BOundaryConditionIdx)%BoundaryConditions, &
                & all_DependentField(BOundaryConditionIdx)%DependentField, &
                & CMFE_FIELD_U_VARIABLE_TYPE,1,1,NodeNumber, ComponentIdx, &
                & MATCH_BC(BoundaryDirichelet(5+6*(i-1),BOundaryConditionIdx)), &
                & BOundaryConditionValue ,Err)
            ELSEIF (TRIM(BoundaryDirichelet(6+6*(i-1),BOundaryConditionIdx)) == "ADD") THEN
              CALL cmfe_BoundaryConditions_AddNode( &
                & all_BoundaryConditions(BOundaryConditionIdx)%BoundaryConditions, &
                & all_DependentField(BOundaryConditionIdx)%DependentField, &
                & CMFE_FIELD_U_VARIABLE_TYPE,1,1,NodeNumber, ComponentIdx, &
                & MATCH_BC(BoundaryDirichelet(5+6*(i-1),BOundaryConditionIdx)), &
                & BOundaryConditionValue ,Err)
            ELSE
              CALL HANDLE_ERROR("Invalid specifier " &
                & //TRIM(BoundaryDirichelet(6+6*(i-1),BOundaryConditionIdx)) &
                & //" in input file. Expected ADD or SET!")
            END IF
            END IF  ! (TRIM(Constraint(ComponentIdx:ComponentIdx)) == "1")
          END DO  ! ComponentIdx=1,3
        END IF ! NodeDOmain==ComputationalNodeNumber
      END DO  !! nodeIdx=1,SIZE(SurfaceNodes,1)

      DEALLOCATE(SurfaceNodes)

    END DO  !! i = 1,NumberOfDirichelet
    CALL cmfe_SolverEquations_BoundaryConditionsCreateFinish(all_SolverEquations(BOundaryConditionIdx)%SolverEquations,Err)
  END DO ! BOundaryConditionIdx
      !!!!!!


  

    
    
  DO ProblemIdx = 1, NumberOfProblem

    !Solve problem
   
    CALL cmfe_Problem_Solve(all_Problem(ProblemIdx)%Problem,Err)

  END DO  ! ProblemIdx

  !Output solution
  DO FieldIdx = 1,NumberOfFields

    CALL cmfe_Fields_Initialise(all_Fields(FieldIdx)%Fields,Err)
    CALL cmfe_Fields_Create(all_Region(FieldIdx)%Region,all_Fields(FieldIdx)%Fields,Err)
    CALL cmfe_Fields_NodesExport(all_Fields(FieldIdx)%Fields,all_Output(FieldIdx)%NodeExport(1) ,"FORTRAN",Err)
    CALL cmfe_Fields_ElementsExport(all_Fields(FieldIdx)%Fields,all_Output(FieldIdx)%ElementExport(1),"FORTRAN",Err)
    CALL cmfe_Fields_Finalise(all_Fields(FieldIdx)%Fields,Err)

  END DO ! FieldIdx

   !!!!!!!!!!!!!!!! IN PROCESS OF IMPLEMENTING IDEA !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!  CALL cmfe_GeneratedMesh_SurfaceGet(all_GeneratedMesh(1)%GeneratedMesh, &
!   & CMFE_GENERATED_MESH_REGULAR_LEFT_SURFACE,LeftSurfaceNodes,BottomNormalXi,Err)

!  CALL cmfe_GeneratedMesh_SurfaceGet(all_GeneratedMesh(1)%GeneratedMesh, &
!    & CMFE_GENERATED_MESH_REGULAR_LEFT_SURFACE,FrontSurfaceNodes,BottomNormalXi,Err)

  
!  !Set z=0 nodes to no z displacementcmgui 
!  DO t=1,SIZE(BottomSurfaceNodes,1)
!    NodeNumber=BottomSurfaceNodes(t)
!  ENDDO
  !!!!!!!!!!!!!!!!!!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  CALL cmfe_Finalise(Err)

  WRITE(*,'(A)') "Program successfully completed."

   

  STOP

  CONTAINS

  INCLUDE "problems_match_library.f90"
  

  !!!!!!!!!!!!!!!! THE FOLLOWING SUBROUTINE IS  USED FOR PRESCRIBING QUADRATIC BOUDARY CONDITIONS !!!!!!!!!!!!!
  SUBROUTINE QUADRATIC(x,y,z, A, B,C,D,E,F,G,H,I,J, BoundaryCOnditionValue)

    REAL(CMISSRP) ,intent(in) :: x,y,z   !!! spatial coordinates
    REAL(CMISSRP) ,intent(in) :: A, B,C,D,E,F,G,H,I,J  !! constants
    REAL(CMISSRP) ,intent(out):: BoundaryCOnditionValue   !!! value of the output functions

    BoundaryCOnditionValue = A*x**2 + B*y**2 + C*z**2 + D*x*y + E*x*z + F*y*z  &
                               & + G*x + H*y + I*z + J


  END SUBROUTINE QUADRATIC

  

END PROGRAM GENERICEXAMPLE
