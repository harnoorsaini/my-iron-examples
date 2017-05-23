! ---------------------------------------------------------------------------------------------------------------------------------
!
! SIMULATION OF A 3D HILL TYPE MUSCLE MODEL
!
! - Displacement controlled - isometric or concentric
! - Muscle activation - variation with time (read in from text file)
! - Muscle activation - variation over space
! - Generated (structed) mesh
!
!
! AUTHOR: HARNOOR SAINI
! MARCH 2017
!
!
! ------------------------ ORIGINAL FILE HEADER START ----------------------------
!> \file
!> \author Adam Reeve
!> \brief This is an example program to solve a finite elasticity equation using OpenCMISS calls.
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
!> Kingdom. Portions created by the University of Auckland and University
!> of Oxford are Copyright (C) 2007 by the University of Auckland and
!> the University of Oxford. All Rights Reserved.
!>
!> Contributor(s): Adam Reeve, Thomas Heidlauf
!>
!> Alternatively, the contents of this file may be used under the terms of
!> either the GNU General Public License Version 2 or later (the "GPL"), or
!> the GNU Lesser General Public License Version 2.1 or later (the "LGPL"),
!> in which case the provisions of the GPL or the LGPL are applicable instead
!> of those above. If you wish to allow use of your version of this file only
!> under the terms of either the GPL or the LGPL, and not to allow others to
!> use your version of this file under the terms of the MPL, indicate your
!> decision by deleting the provisions above and replace them with the notice
!> and other provisions required by the GPL or the LGPL. If you do not delete
!> the provisions above, a recipient may use your version of this file under
!> the terms of any one of the MPL, the GPL or the LGPL.
! ------------------------ ORIGINAL FILE HEADER END ---------------------------- 
!

!> Main program
PROGRAM LARGEUNIAXIALEXTENSIONEXAMPLE

  USE OpenCMISS
  USE OpenCMISS_Iron
#ifndef NOMPIMOD
  USE MPI
#endif

#ifdef WIN32
  USE IFQWIN
#endif

  IMPLICIT NONE

#ifdef NOMPIMOD
#include "mpif.h"
#endif

!
! Parameter Definition--------------------------------------------------------------------------------------------------------START
!

  ! Equations set for transversely isotropic (fibre-reinforced), active contractible bodies consitisting of two materials
  ! The local portion between them is defined by the parameter trans
  ! Material 1 is active contractible, material 2 is only passive
  !   W=W_iso+W_aniso+W_act
  !     where the three parts are adopted from above (iso Mooney-Rivlin, aniso Markert, active part)
  !   Markert, B., W. Ehlers, and N. Karajan. 
  !   A general polyconvex strain-energy function for fiber-reinforced materials. 
  !   Proceedings in Applied Mathematics and Mechanics 5.1 (2005): 245-246.)
  ! 
  ! C(1)=c1_m1...Mooney Rivlin parameter material 1
  ! C(2)=c2_m1...Mooney Rivlin parameter material 1
  ! C(3)=c4_m1...polynomial coefficient (Markert model) material 1
  ! C(4)=c5_m1...power coefficient (Markert model) material 1
  ! C(5)=c1_m2...Mooney Rivlin parameter material 2
  ! C(6)=c2_m2...Mooney Rivlin parameter material 2
  ! C(7)=c4_m2...polynomial coefficient (Markert model) material 2
  ! C(8)=c5_m2...power coefficient (Markert model) material 2
  ! C(9)=alpha...activation parameter [0,1]
  ! C(10)=trans...transition parameter [0,1] for the portion between the two materials
  ! C(11)=P_max...maximum isometric stress

  REAL(CMISSRP), PARAMETER, DIMENSION(11) :: C= &
    & [3.56E-2_CMISSRP,3.86E-2_CMISSRP,0.3E-8_CMISSRP, &
    &  34.0_CMISSRP,3.56E-2_CMISSRP, 3.86E-2_CMISSRP, 0.3E-8_CMISSRP, &
    &  34.0_CMISSRP, 0.0_CMISSRP, 1.0_CMISSRP, 0.15_CMISSRP] 

    ! From original Biceps example file
    !REAL(CMISSRP), PARAMETER, DIMENSION(11) :: C= &
    !& [0.0356_CMISSRP,0.00386_CMISSRP, &
    !& 0.0000000357_CMISSRP,42.6_CMISSRP, &
    !& 2.31_CMISSRP,0.00000115_CMISSRP, &
    !& 7.99_CMISSRP,16.6_CMISSRP, &
    !& 0.0_CMISSRP,1.0_CMISSRP, &
    !& 0.3_CMISSRP] 

  REAL(CMISSRP), PARAMETER :: gamma=0.01_CMISSRP, beta=1.0_CMISSRP
    
  ! Test program parameters
  REAL(CMISSDP), PARAMETER :: PI=4.0_CMISSDP*DATAN(1.0_CMISSDP)

  REAL(CMISSRP), PARAMETER :: HEIGHT=2.0_CMISSRP
  REAL(CMISSRP), PARAMETER :: WIDTH=5.0_CMISSRP
  REAL(CMISSRP), PARAMETER :: LENGTH=1.0_CMISSRP
  !INTEGER(CMISSIntg), PARAMETER :: InterpolationType=CMFE_BASIS_LINEAR_LAGRANGE_INTERPOLATION
  INTEGER(CMISSIntg), PARAMETER :: InterpolationType=CMFE_BASIS_QUADRATIC_LAGRANGE_INTERPOLATION
  INTEGER(CMISSIntg), PARAMETER :: PressureInterpolationType=CMFE_BASIS_LINEAR_LAGRANGE_INTERPOLATION
  LOGICAL, PARAMETER :: UsePressureBasis=.TRUE.
  INTEGER(CMISSIntg), PARAMETER :: NumberOfGaussXi=3

  INTEGER(CMISSIntg), PARAMETER :: CoordinateSystemUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: RegionUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: BasisUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: PressureBasisUserNumber=2
  INTEGER(CMISSIntg), PARAMETER :: GeneratedMeshUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: MeshUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: DecompositionUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: FieldGeometryUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: FieldFibreUserNumber=2
  INTEGER(CMISSIntg), PARAMETER :: FieldMaterialUserNumber=3
  INTEGER(CMISSIntg), PARAMETER :: FieldDependentUserNumber=4
  INTEGER(CMISSIntg), PARAMETER :: EquationSetUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: EquationsSetFieldUserNumber=5
  INTEGER(CMISSIntg), PARAMETER :: ProblemUserNumber=1
  
  INTEGER(CMISSIntg), PARAMETER :: NumberOfElements=39
  INTEGER(CMISSIntg), PARAMETER :: NumberOfNodes=567
  INTEGER(CMISSIntg), PARAMETER :: OutputFrequency=1

  INTEGER(CMISSIntg), DIMENSION(8), PARAMETER :: Entries=[1,3,7,9,19,21,25,27]

  INTEGER(CMISSIntg), DIMENSION(21), PARAMETER :: TendonTopFixNodes=[35,36,37,38,41,42,43,44,119,120,121,122,139,140,179,180, &
    & 181,182,199,200,205]
  INTEGER(CMISSIntg), DIMENSION(21), PARAMETER :: TendonBottomFixNodes=[5,6,7,8,21,22,73,74,89,90,91,92,127,128,149,150,151, &
    & 152,187,188,208]


  ! Program types
  ! Program variables
  INTEGER(CMISSIntg) :: NumberGlobalXElements,NumberGlobalYElements,NumberGlobalZElements
  INTEGER(CMISSIntg) :: EquationsSetIndex, TotalNumberOfNodes
  INTEGER(CMISSIntg) :: NumberOfComputationalNodes,NumberOfDomains,ComputationalNodeNumber
  INTEGER(CMISSIntg) :: NodeNumber,NodeDomain,node_idx,elem_idx, gauss_idx, component_idx, source_idx

  INTEGER(CMISSIntg) :: FT_1,FT_2,FT_3,FT_4,FT_5,FT_6,FT_7,FT_8,FT_9, &
    & FT_10,FT_11,FT_12,FT_13,FT_14,FT_15,FT_16,FT_17,FT_18,FT_19, &
    & FT_20,FT_21,FT_22,FT_23,FT_24,FT_25,FT_26,FT_27
  
  REAL(CMISSRP) :: x1,x2,x3
  REAL(CMISSRP) :: VALUE,INIT_PRESSURE

  REAL(CMISSRP), DIMENSION(NumberOfNodes,3) :: AllNodes
  INTEGER(CMISSIntg), DIMENSION(NumberOfElements,27) :: AllElements

  INTEGER(CMISSIntg) :: Elem,Node
  INTEGER(CMISSIntg) :: stat
  character(len=256) :: filename,string
  ! to delete
  INTEGER(CMISSIntg),ALLOCATABLE :: BottomSurfaceNodes(:)
  INTEGER(CMISSIntg),ALLOCATABLE :: LeftSurfaceNodes(:)
  INTEGER(CMISSIntg),ALLOCATABLE :: RightSurfaceNodes(:)
  INTEGER(CMISSIntg),ALLOCATABLE :: FrontSurfaceNodes(:)
  INTEGER(CMISSIntg),ALLOCATABLE :: TopSurfaceNodes(:)
  INTEGER(CMISSIntg),ALLOCATABLE :: BackSurfaceNodes(:)
  ! to delete
  INTEGER(CMISSIntg) :: BottomNormalXi,LeftNormalXi,RightNormalXi,FrontNormalXi
  INTEGER(CMISSIntg) :: NumberOfElementsFE, FibreFieldNumberOfComponents
  !INTEGER(CMISSIntg), PARAMETER :: NUMBER_OF_COMPONENTS = 3 !nearly incompressible
  INTEGER(CMISSIntg), PARAMETER :: NUMBER_OF_COMPONENTS = 4 !fully incompressible

  ! CMISS variables
  TYPE(cmfe_BasisType) :: Basis, PressureBasis
  TYPE(cmfe_BoundaryConditionsType) :: BoundaryConditions
  TYPE(cmfe_CoordinateSystemType) :: CoordinateSystem, WorldCoordinateSystem
  TYPE(cmfe_MeshType) :: Mesh
  TYPE(cmfe_DecompositionType) :: Decomposition
  TYPE(cmfe_EquationsType) :: Equations
  TYPE(cmfe_EquationsSetType) :: EquationsSet
  TYPE(cmfe_FieldType) :: GeometricField,MaterialField,DependentField,EquationsSetField,FibreField
  TYPE(cmfe_FieldsType) :: Fields
  TYPE(cmfe_ProblemType) :: Problem
  TYPE(cmfe_RegionType) :: Region,WorldRegion
  TYPE(cmfe_SolverType) :: Solver,MechanicalSolver
  TYPE(cmfe_SolverEquationsType) :: SolverEquations
  TYPE(cmfe_ControlLoopType) :: ControlLoop
  TYPE(cmfe_NodesType) :: Nodes
  TYPE(cmfe_MeshElementsType) :: QuadraticElements
  TYPE(cmfe_MeshElementsType) :: LinearElements

  INTEGER(CMISSIntg) :: i, k, j
  
#ifdef WIN32
  !Quickwin type
  LOGICAL :: QUICKWIN_STATUS=.FALSE.
  TYPE(WINDOWCONFIG) :: QUICKWIN_WINDOW_CONFIG
#endif

  ! Generic CMISS variables
  INTEGER(CMISSIntg) :: Err

  ! Variables, Parameters, ...for specific simulation
  INTEGER(CMISSIntg), PARAMETER :: TIMESTEPS=10 ! Number of Timesteps
  INTEGER(CMISSIntg), PARAMETER :: TotalNumberOfSources=2
  INTEGER(CMISSIntg) :: alpha_SU ! spatial update of alpha
  REAL(CMISSRP), DIMENSION(TIMESTEPS) :: alpha
  REAL(CMISSRP) :: alpha_t, BCLOAD, alpha_spatial_temporal, alpha_spatial_wt_sum
  REAL(CMISSRP) :: FibreFieldAngle(3), nodal_Coords(3)
  REAL(CMISSRP), ALLOCATABLE :: dist_Node_actSource(:,:), alpha_spatial_wt(:,:)
  REAL(CMISSRP), DIMENSION(3,TotalNumberOfSources) :: activation_Source ! NOTE: dimensions will be transposed by read command
!
! Initialize parameters and set up mesh---------------------------------------------------------------------------------------START
!

#ifdef WIN32
  ! Initialise QuickWin
  QUICKWIN_WINDOW_CONFIG%TITLE="General Output" !Window title
  QUICKWIN_WINDOW_CONFIG%NUMTEXTROWS=-1 !Max possible number of rows
  QUICKWIN_WINDOW_CONFIG%MODE=QWIN$SCROLLDOWN
  ! Set the window parameters
  QUICKWIN_STATUS=SETWINDOWCONFIG(QUICKWIN_WINDOW_CONFIG)
  ! If attempt fails set with system estimated values
  IF(.NOT.QUICKWIN_STATUS) QUICKWIN_STATUS=SETWINDOWCONFIG(QUICKWIN_WINDOW_CONFIG)
#endif

  ! Read in mesh ---------------------------------------------------------------------------------------START
 WRITE(*,*) "Reading file: input/biceps_3D27.dat"
  OPEN(UNIT=1,FILE="input/biceps_3D27_new.dat",IOSTAT=stat)

  !read the coodinates of the nodes from file
  k=1
  READ(1,*,IOSTAT=stat) 
  READ(1,*,IOSTAT=stat)
  DO WHILE(.true.)
    READ(1,*,IOSTAT=stat) Node,x1,x2,x3
    IF(stat<0) EXIT !end of file
    IF(k/=Node) THEN
      WRITE(*,*) "Error in reading nodal coordinates"
    END IF
    
    AllNodes(k,:)=[x1,x2,x3]
    IF(k==NumberOfNodes) THEN
      WRITE(*,*) "There are ", NumberOfNodes, " nodes in the mesh"
      EXIT
    END IF
    k=k+1
  END DO

  !read the element nodes from file
  k=1
  READ(1,*,IOSTAT=stat) 
  READ(1,*,IOSTAT=stat) 
  READ(1,*,IOSTAT=stat)
  DO WHILE(.true.)
    READ(1,*,IOSTAT=stat) Elem,FT_1,FT_2,FT_3,FT_4,FT_5,FT_6,FT_7,FT_8,FT_9, &
      & FT_10,FT_11,FT_12,FT_13,FT_14,FT_15,FT_16,FT_17,FT_18,FT_19, &
      & FT_20,FT_21,FT_22,FT_23,FT_24,FT_25,FT_26,FT_27

    IF(stat<0) EXIT !end of file
    IF(k/=Elem) THEN
      WRITE(*,*) "Error in reading element nodes"
      EXIT
    END IF
    
    !change order such that Xi_1 rather than Xi_3 direction is aligned with the fibre direction
    AllElements(k,:)=[FT_19,FT_10,FT_1,FT_22,FT_13,FT_4,FT_25,FT_16,FT_7, &
      & FT_20,FT_11,FT_2,FT_23,FT_14,FT_5,FT_26,FT_17,FT_8, &
      & FT_21,FT_12,FT_3,FT_24,FT_15,FT_6,FT_27,FT_18,FT_9]
    IF(k==NumberOfElements) THEN
      WRITE(*,*) "There are ", NumberOfElements, " elements in the mesh"
      EXIT
    END IF
    k=k+1
  END DO

  CLOSE(UNIT=1)
  WRITE(*,*) "Finished reading file: input/biceps_3D27.dat"
  ! Read in mesh -----------------------------------------------------------------------------------------END

  ! Intialise cmiss
  CALL cmfe_Initialise(WorldCoordinateSystem,WorldRegion,Err)
  CALL cmfe_ErrorHandlingModeSet(CMFE_ERRORS_TRAP_ERROR,Err)

  ! Set all diganostic levels on for testing
  !CALL cmfe_DiagnosticsSetOn(CMFE_FROM_DIAG_TYPE,[1,2,3,4,5],"Diagnostics",["PROBLEM_RESIDUAL_EVALUATE"],Err)

  ! Get the number of computational nodes and this computational node number
  CALL cmfe_ComputationalNumberOfNodesGet(NumberOfComputationalNodes,Err)
  CALL cmfe_ComputationalNodeNumberGet(ComputationalNodeNumber,Err)
  
  WRITE(*,*) "Comp-Node Number", "-", "Total Number of Comp-Nodes"  
  WRITE(*,*) ComputationalNodeNumber, NumberOfComputationalNodes

  !NumberGlobalXElements=5
  !NumberGlobalYElements=4
  !NumberGlobalZElements=2
  !NumberOfElementsFE=NumberGlobalXElements*NumberGlobalYElements*NumberGlobalZElements
  NumberOfDomains=NumberOfComputationalNodes

  ! Create a 3D rectangular cartesian coordinate system
  CALL cmfe_CoordinateSystem_Initialise(CoordinateSystem,Err)
  CALL cmfe_CoordinateSystem_CreateStart(CoordinateSystemUserNumber,CoordinateSystem,Err)
  CALL cmfe_CoordinateSystem_CreateFinish(CoordinateSystem,Err)

  ! Create a region and assign the coordinate system to the region
  CALL cmfe_Region_Initialise(Region,Err)
  CALL cmfe_Region_CreateStart(RegionUserNumber,WorldRegion,Region,Err)
  CALL cmfe_Region_LabelSet(Region,"Region",Err)
  CALL cmfe_Region_CoordinateSystemSet(Region,CoordinateSystem,Err)
  CALL cmfe_Region_CreateFinish(Region,Err)

!
! Basis-----------------------------------------------------------------------------------------------------------------------START
!

  !Define geometric basis
  CALL cmfe_Basis_Initialise(Basis,Err)
  CALL cmfe_Basis_CreateStart(BasisUserNumber,Basis,Err)
  CALL cmfe_Basis_TypeSet(Basis,CMFE_BASIS_LAGRANGE_HERMITE_TP_TYPE,Err)
  CALL cmfe_Basis_NumberOfXiSet(Basis,3,Err)
  CALL cmfe_Basis_InterpolationXiSet(Basis,[CMFE_BASIS_QUADRATIC_LAGRANGE_INTERPOLATION, &
    & CMFE_BASIS_QUADRATIC_LAGRANGE_INTERPOLATION,CMFE_BASIS_QUADRATIC_LAGRANGE_INTERPOLATION],Err)
  CALL cmfe_Basis_QuadratureNumberOfGaussXiSet(Basis,[NumberOfGaussXi,NumberOfGaussXi,NumberOfGaussXi],Err)
  CALL cmfe_Basis_CreateFinish(Basis,Err)

  !Define pressure basis
  CALL cmfe_Basis_Initialise(PressureBasis,Err)
  CALL cmfe_Basis_CreateStart(PressureBasisUserNumber,PressureBasis,Err)
  CALL cmfe_Basis_TypeSet(PressureBasis,CMFE_BASIS_LAGRANGE_HERMITE_TP_TYPE,Err)
  CALL cmfe_Basis_NumberOfXiSet(PressureBasis,3,Err)
  CALL cmfe_Basis_InterpolationXiSet(PressureBasis,[CMFE_BASIS_LINEAR_LAGRANGE_INTERPOLATION, &
    & CMFE_BASIS_LINEAR_LAGRANGE_INTERPOLATION,CMFE_BASIS_LINEAR_LAGRANGE_INTERPOLATION],Err)
  CALL cmfe_Basis_QuadratureNumberOfGaussXiSet(PressureBasis,[NumberOfGaussXi,NumberOfGaussXi,NumberOfGaussXi],Err)
  CALL cmfe_Basis_CreateFinish(PressureBasis,Err)

!
! Meshing---------------------------------------------------------------------------------------------------------------------START
!
  !Create a mesh with three-dimensional elements
  CALL cmfe_Mesh_Initialise(Mesh,Err)
  CALL cmfe_Mesh_CreateStart(MeshUserNumber,Region,3,Mesh,Err)
  CALL cmfe_Mesh_NumberOfComponentsSet(Mesh,2,Err) 
  CALL cmfe_Mesh_NumberOfElementsSet(Mesh,NumberOfElements,Err)  
  !Define nodes for the mesh
  CALL cmfe_Nodes_Initialise(Nodes,Err)
  CALL cmfe_Nodes_CreateStart(Region,NumberOfNodes,Nodes,Err)
  CALL cmfe_Nodes_CreateFinish(Nodes,Err)

  CALL cmfe_MeshElements_Initialise(QuadraticElements,Err)
  CALL cmfe_MeshElements_CreateStart(Mesh,1,Basis,QuadraticElements,Err)
  DO j=1,NumberOfElements
    CALL cmfe_MeshElements_NodesSet(QuadraticElements,j,AllElements(j,:),Err)
  ENDDO
  CALL cmfe_MeshElements_CreateFinish(QuadraticElements,Err)

  !for the linear elements, we need only specific entries from the elements above
  CALL cmfe_MeshElements_Initialise(LinearElements,Err)
  CALL cmfe_MeshElements_CreateStart(Mesh,2,PressureBasis,LinearElements,Err)
  DO j=1,NumberOfElements
    CALL cmfe_MeshElements_NodesSet(LinearElements,j,AllElements(j,Entries),Err)
  ENDDO
  CALL cmfe_MeshElements_CreateFinish(LinearElements,Err)

  CALL cmfe_Mesh_CreateFinish(Mesh,Err) 
  ! Start the creation of a generated mesh in the region
  !CALL cmfe_GeneratedMesh_Initialise(GeneratedMesh,Err)
  !CALL cmfe_GeneratedMesh_CreateStart(GeneratedMeshUserNumber,Region,GeneratedMesh,Err)
  ! Set up a regular x*y*z mesh
  !CALL cmfe_GeneratedMesh_TypeSet(GeneratedMesh,CMFE_GENERATED_MESH_REGULAR_MESH_TYPE,Err)
  ! Set the default basis
  !IF(UsePressureBasis) THEN
    !CALL cmfe_GeneratedMesh_BasisSet(GeneratedMesh,[Basis,PressureBasis],Err)
  !ELSE
    !CALL cmfe_GeneratedMesh_BasisSet(GeneratedMesh,[Basis],Err)
  !ENDIF
  ! Define the mesh on the region
  !IF(NumberGlobalXElements==0) THEN
    !CALL cmfe_GeneratedMesh_ExtentSet(GeneratedMesh,[WIDTH,HEIGHT],Err)
    !CALL cmfe_GeneratedMesh_NumberOfElementsSet(GeneratedMesh,[NumberGlobalXElements,NumberGlobalYElements],Err)
  !ELSE
    !CALL cmfe_GeneratedMesh_ExtentSet(GeneratedMesh,[WIDTH,HEIGHT,LENGTH],Err)
    !CALL cmfe_GeneratedMesh_NumberOfElementsSet(GeneratedMesh,[NumberGlobalXElements,NumberGlobalYElements, &
     !& NumberGlobalZElements],Err)
  !ENDIF
  ! Finish the creation of a generated mesh in the region
  !CALL cmfe_Mesh_Initialise(Mesh,Err)
  !CALL cmfe_GeneratedMesh_CreateFinish(GeneratedMesh,MeshUserNumber,Mesh,Err)

  ! Create a decomposition
  CALL cmfe_Decomposition_Initialise(Decomposition,Err)
  CALL cmfe_Decomposition_CreateStart(DecompositionUserNumber,Mesh,Decomposition,Err)
  CALL cmfe_Decomposition_TypeSet(Decomposition,CMFE_DECOMPOSITION_CALCULATED_TYPE,Err)
  CALL cmfe_Decomposition_NumberOfDomainsSet(Decomposition,NumberOfDomains,Err)
  CALL cmfe_Decomposition_CreateFinish(Decomposition,Err)

!
! Define fields-------------------------------------------------------------------------------------------------------------START
!

  ! Create a field to put the geometry (default is geometry)
  CALL cmfe_Field_Initialise(GeometricField,Err)
  CALL cmfe_Field_CreateStart(FieldGeometryUserNumber,Region,GeometricField,Err)
  CALL cmfe_Field_MeshDecompositionSet(GeometricField,Decomposition,Err)
  CALL cmfe_Field_VariableLabelSet(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,"Geometry",Err)
  CALL cmfe_Field_ScalingTypeSet(GeometricField,CMFE_FIELD_ARITHMETIC_MEAN_SCALING,Err)
  CALL cmfe_Field_CreateFinish(GeometricField,Err)

  !initialise the Geometry
  DO node_idx=1,SIZE(AllNodes,1)  
    CALL cmfe_Decomposition_NodeDomainGet(Decomposition,node_idx,1,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN
      CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1, &
       & node_idx,1,AllNodes(node_idx,1),Err)
      CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1, &
       & node_idx,2,AllNodes(node_idx,2),Err)
      CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1, &
       & node_idx,3,AllNodes(node_idx,3),Err)
    END IF
  END DO
  ! Update the geometric field parameters
  !CALL cmfe_GeneratedMesh_GeometricParametersCalculate(GeneratedMesh,GeometricField,Err)

  ! Create a fibre field and attach it to the geometric field
  CALL cmfe_Field_Initialise(FibreField,Err)
  CALL cmfe_Field_CreateStart(FieldFibreUserNumber,Region,FibreField,Err)
  CALL cmfe_Field_TypeSet(FibreField,CMFE_FIELD_FIBRE_TYPE,Err)
  CALL cmfe_Field_MeshDecompositionSet(FibreField,Decomposition,Err)
  CALL cmfe_Field_GeometricFieldSet(FibreField,GeometricField,Err)
  CALL cmfe_Field_VariableLabelSet(FibreField,CMFE_FIELD_U_VARIABLE_TYPE,"Fibre",Err)

  ! >>>> ADDED IN FROM MYLENA FOR FIBER ORIENTATION >>>> START
  CALL cmfe_Field_NumberOfVariablesSet(FibreField,1,Err)
  FibreFieldNumberOfComponents=3
  CALL cmfe_Field_NumberOfComponentsSet(FibreField,CMFE_FIELD_U_VARIABLE_TYPE,FibreFieldNumberOfComponents,Err) 
  CALL cmfe_Field_ComponentMeshComponentSet(FibreField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(FibreField,CMFE_FIELD_U_VARIABLE_TYPE,2,1,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(FibreField,CMFE_FIELD_U_VARIABLE_TYPE,3,1,Err)
  
  !                                                 >>>> END
  CALL cmfe_Nodes_Initialise(Nodes,Err)
  CALL cmfe_Region_NodesGet(Region,Nodes,Err)
  CALL cmfe_Nodes_NumberOfNodesGet(Nodes,TotalNumberOfNodes,Err)
  CALL cmfe_Field_CreateFinish(FibreField,Err)

  ! >>>> ADDED IN FROM MYLENA FOR FIBER ORIENTATION >>>> START
  ! Rotation Angles (in radiant!!)
  ! in 2D an entry in Angle(1) means rotated x-axis, 
  !          entry in Angle(2) doesn't make sense, as rotates out of surface ...
  ! in 3D an entry in Angle(1) means rotated around z-axis, entry in Angle(2) means rotated around y-axis
  !          entry in Angle(3) means rotated around x-axis => no change
  ! 45° equivalent to pi/4, 90° equivalent to pi/2
  ! 1 degree = 0.0174533rad

  ! FibreFieldAngle=(/10.0_CMISSRP*0.0174533_CMISSRP,0.0_CMISSRP,0.0_CMISSRP/)
  FibreFieldAngle=(/0.0_CMISSRP,0.0_CMISSRP,0.0_CMISSRP/)

  DO node_idx=1,TotalNumberOfNodes
    CALL cmfe_Decomposition_NodeDomainGet(Decomposition,node_idx,1,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN
      DO component_idx=1,FibreFieldNumberOfComponents
        CALL cmfe_Field_ParameterSetUpdateNode(FibreField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1, &
          & 1,node_idx,component_idx,FibreFieldAngle(component_idx),Err)
      ENDDO
    ENDIF
  ENDDO
  !                                                 >>>> END

  ! Create the material field
  CALL cmfe_Field_Initialise(MaterialField,Err)
  CALL cmfe_Field_CreateStart(FieldMaterialUserNumber,Region,MaterialField,Err)
  CALL cmfe_Field_TypeSet(MaterialField,CMFE_FIELD_MATERIAL_TYPE,Err)
  CALL cmfe_Field_MeshDecompositionSet(MaterialField,Decomposition,Err)
  CALL cmfe_Field_GeometricFieldSet(MaterialField,GeometricField,Err)
  CALL cmfe_Field_NumberOfVariablesSet(MaterialField,1,Err)
  CALL cmfe_Field_VariableLabelSet(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,"Material",Err)
  CALL cmfe_Field_NumberOfComponentsSet(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,11,Err)
  CALL cmfe_Field_ComponentInterpolationSet(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,1,CMFE_FIELD_ELEMENT_BASED_INTERPOLATION,Err)
  CALL cmfe_Field_ComponentInterpolationSet(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,2,CMFE_FIELD_ELEMENT_BASED_INTERPOLATION,Err)
  CALL cmfe_Field_ComponentInterpolationSet(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,3,CMFE_FIELD_ELEMENT_BASED_INTERPOLATION,Err)
  CALL cmfe_Field_ComponentInterpolationSet(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,4,CMFE_FIELD_ELEMENT_BASED_INTERPOLATION,Err)
  ! Spatial variation of muscle activation 
    ! Required since value of activation is different at each point (see line 3303 in opencmiss_iron.f90)
    
    ! Element based interpolation
    !CALL cmfe_Field_ComponentInterpolationSet(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,5,CMFE_FIELD_ELEMENT_BASED_INTERPOLATION,Err)
    
    ! Node based interpolation
    CALL cmfe_Field_ComponentInterpolationSet(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,5,CMFE_FIELD_ELEMENT_BASED_INTERPOLATION,Err)
    
    ! Gauss-point based interpolation
    !CALL cmfe_Field_ComponentInterpolationSet(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,5, & 
    !& CMFE_FIELD_GAUSS_POINT_BASED_INTERPOLATION,Err)
  CALL cmfe_Field_ComponentInterpolationSet(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,6,CMFE_FIELD_ELEMENT_BASED_INTERPOLATION,Err)
  CALL cmfe_Field_ComponentInterpolationSet(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,7,CMFE_FIELD_ELEMENT_BASED_INTERPOLATION,Err)
  CALL cmfe_Field_ComponentInterpolationSet(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,8,CMFE_FIELD_ELEMENT_BASED_INTERPOLATION,Err)
  CALL cmfe_Field_ComponentInterpolationSet(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,9,CMFE_FIELD_ELEMENT_BASED_INTERPOLATION,Err)
  CALL cmfe_Field_ComponentInterpolationSet(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,10,CMFE_FIELD_ELEMENT_BASED_INTERPOLATION,Err)
  CALL cmfe_Field_ComponentInterpolationSet(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,11,CMFE_FIELD_ELEMENT_BASED_INTERPOLATION,Err)
  CALL cmfe_Field_CreateFinish(MaterialField,Err)

  ! Create the dependent field
  CALL cmfe_Field_Initialise(DependentField,Err)
  CALL cmfe_Field_CreateStart(FieldDependentUserNumber,Region,DependentField,Err)
  CALL cmfe_Field_TypeSet(DependentField,CMFE_FIELD_GEOMETRIC_GENERAL_TYPE,Err)
  CALL cmfe_Field_MeshDecompositionSet(DependentField,Decomposition,Err)
  CALL cmfe_Field_GeometricFieldSet(DependentField,GeometricField,Err)
  CALL cmfe_Field_DependentTypeSet(DependentField,CMFE_FIELD_DEPENDENT_TYPE,Err)
  CALL cmfe_Field_NumberOfVariablesSet(DependentField,2,Err)
  CALL cmfe_Field_VariableTypesSet(DependentField,[CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_DELUDELN_VARIABLE_TYPE],Err)
  CALL cmfe_Field_VariableLabelSet(DependentField,CMFE_FIELD_U_VARIABLE_TYPE,"Dependent",Err)
  CALL cmfe_Field_NumberOfComponentsSet(DependentField,CMFE_FIELD_U_VARIABLE_TYPE,NUMBER_OF_COMPONENTS,Err)
  CALL cmfe_Field_NumberOfComponentsSet(DependentField,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,NUMBER_OF_COMPONENTS,Err)
  !CALL cmfe_Field_NumberOfComponentsSet(DependentField,CMFE_FIELD_V_VARIABLE_TYPE,3,Err)
  IF(UsePressureBasis) THEN
    ! Set the pressure to be nodally based and use the second mesh component if required
    CALL cmfe_Field_ComponentInterpolationSet(DependentField,CMFE_FIELD_U_VARIABLE_TYPE,4,CMFE_FIELD_NODE_BASED_INTERPOLATION,Err)
    CALL cmfe_Field_ComponentInterpolationSet(DependentField,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,4, &
      & CMFE_FIELD_NODE_BASED_INTERPOLATION,Err)
    CALL cmfe_Field_ComponentMeshComponentSet(DependentField,CMFE_FIELD_U_VARIABLE_TYPE,4,2,Err)
    CALL cmfe_Field_ComponentMeshComponentSet(DependentField,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,4,2,Err)
  END IF
  CALL cmfe_Field_CreateFinish(DependentField,Err)

!
! Equaton set-----------------------------------------------------------------------------------------------------------------START
!

  ! Create the equations_set
  CALL cmfe_Field_Initialise(EquationsSetField,Err)
  CALL cmfe_EquationsSet_CreateStart(EquationSetUserNumber,Region,FibreField,[CMFE_EQUATIONS_SET_ELASTICITY_CLASS, &
    & CMFE_EQUATIONS_SET_FINITE_ELASTICITY_TYPE,CMFE_EQUATIONS_SET_TRANS_ISOTROPIC_ACTIVE_TRANSITION_SUBTYPE], &
    & EquationsSetFieldUserNumber,EquationsSetField,EquationsSet,Err)
  CALL cmfe_EquationsSet_CreateFinish(EquationsSet,Err)

  ! Create the equations set dependent field
  CALL cmfe_EquationsSet_DependentCreateStart(EquationsSet,FieldDependentUserNumber,DependentField,Err)
  CALL cmfe_EquationsSet_DependentCreateFinish(EquationsSet,Err)

  ! Create the equations set material field 
  CALL cmfe_EquationsSet_MaterialsCreateStart(EquationsSet,FieldMaterialUserNumber,MaterialField,Err)
  CALL cmfe_EquationsSet_MaterialsCreateFinish(EquationsSet,Err)

  ! Set Material-Parameters [mu(1) mu(2) mu(3) alpha(1) alpha(2) alpha(3) mu_0 XB]
  CALL cmfe_Field_ComponentValuesInitialise(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,C(1),Err)
  CALL cmfe_Field_ComponentValuesInitialise(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,2,C(2),Err)
  CALL cmfe_Field_ComponentValuesInitialise(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,3,C(3),Err)
  CALL cmfe_Field_ComponentValuesInitialise(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,4,C(4),Err)
  CALL cmfe_Field_ComponentValuesInitialise(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,5,C(5),Err)
  CALL cmfe_Field_ComponentValuesInitialise(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,6,C(6),Err)
  CALL cmfe_Field_ComponentValuesInitialise(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,7,C(7),Err)
  CALL cmfe_Field_ComponentValuesInitialise(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,8,C(8),Err)
  CALL cmfe_Field_ComponentValuesInitialise(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,9,C(9),Err)
  CALL cmfe_Field_ComponentValuesInitialise(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,10,C(10),Err)
  CALL cmfe_Field_ComponentValuesInitialise(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,11,C(11),Err)

! fem group elem  5..8,21..22,28..29,32..36,38..39      --> 100% muscle
  ! default -- do nothing

  ! fem group elem  3..4,13..18,20                        --> 50% muscle
  CALL cmfe_Field_ParameterSetUpdateElement(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 3,10,0.5_CMISSRP,Err)
  CALL cmfe_Field_ParameterSetUpdateElement(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 4,10,0.5_CMISSRP,Err)
  CALL cmfe_Field_ParameterSetUpdateElement(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 13,10,0.5_CMISSRP,Err)
  CALL cmfe_Field_ParameterSetUpdateElement(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 14,10,0.5_CMISSRP,Err)
  CALL cmfe_Field_ParameterSetUpdateElement(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 15,10,0.5_CMISSRP,Err)
  CALL cmfe_Field_ParameterSetUpdateElement(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 16,10,0.5_CMISSRP,Err)
  CALL cmfe_Field_ParameterSetUpdateElement(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 17,10,0.5_CMISSRP,Err)
  CALL cmfe_Field_ParameterSetUpdateElement(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 18,10,0.5_CMISSRP,Err)
  CALL cmfe_Field_ParameterSetUpdateElement(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 20,10,0.5_CMISSRP,Err)

  !INIT_PRESSURE=0.5_CMISSRP*(-2.0_CMISSRP*MAT_FE(2)-MAT_FE(1))+0.5_CMISSRP*(-2.0_CMISSRP*MAT_FE(6)-MAT_FE(5))
  !CALL cmfe_Field_ParameterSetUpdateNode(DependentField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1, &
    !& 3,4,INIT_PRESSURE,Err)

  ! fem group elem  11..12                                --> 80% muscle
  CALL cmfe_Field_ParameterSetUpdateElement(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 11,10,0.8_CMISSRP,Err)
  CALL cmfe_Field_ParameterSetUpdateElement(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 12,10,0.8_CMISSRP,Err)

  ! fem group elem  30..31,37                             --> 99% muscle
  CALL cmfe_Field_ParameterSetUpdateElement(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 30,10,0.99_CMISSRP,Err)
  CALL cmfe_Field_ParameterSetUpdateElement(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 31,10,0.99_CMISSRP,Err)
  CALL cmfe_Field_ParameterSetUpdateElement(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 37,10,0.99_CMISSRP,Err)

  ! fem group elem  1..2,9..10,19,23                      --> 99.9% muscle
  CALL cmfe_Field_ParameterSetUpdateElement(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 1,10,0.999_CMISSRP,Err)
  CALL cmfe_Field_ParameterSetUpdateElement(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 2,10,0.999_CMISSRP,Err)
  CALL cmfe_Field_ParameterSetUpdateElement(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 9,10,0.999_CMISSRP,Err)
  CALL cmfe_Field_ParameterSetUpdateElement(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 10,10,0.999_CMISSRP,Err)
  CALL cmfe_Field_ParameterSetUpdateElement(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 19,10,0.999_CMISSRP,Err)
  CALL cmfe_Field_ParameterSetUpdateElement(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 23,10,0.999_CMISSRP,Err)

  ! fem group elem  24..27                                --> 100% soft tissue (all anisotropic contributions = 0)
  CALL cmfe_Field_ParameterSetUpdateElement(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 24,3,0.0_CMISSRP,Err)
  CALL cmfe_Field_ParameterSetUpdateElement(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 24,7,0.0_CMISSRP,Err)
  CALL cmfe_Field_ParameterSetUpdateElement(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 25,3,0.0_CMISSRP,Err)
  CALL cmfe_Field_ParameterSetUpdateElement(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 25,7,0.0_CMISSRP,Err)
  CALL cmfe_Field_ParameterSetUpdateElement(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 26,3,0.0_CMISSRP,Err)
  CALL cmfe_Field_ParameterSetUpdateElement(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 26,7,0.0_CMISSRP,Err)
  CALL cmfe_Field_ParameterSetUpdateElement(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 27,3,0.0_CMISSRP,Err)
  CALL cmfe_Field_ParameterSetUpdateElement(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 27,7,0.0_CMISSRP,Err)

  ! Create the equations set equations
  CALL cmfe_Equations_Initialise(Equations,Err)
  CALL cmfe_EquationsSet_EquationsCreateStart(EquationsSet,Equations,Err)
  CALL cmfe_Equations_SparsityTypeSet(Equations,CMFE_EQUATIONS_SPARSE_MATRICES,Err)
  CALL cmfe_Equations_OutputTypeSet(Equations,CMFE_EQUATIONS_NO_OUTPUT,Err)
  CALL cmfe_EquationsSet_EquationsCreateFinish(EquationsSet,Err)

  ! Initialise dependent field from undeformed geometry and displacement bcs and set hydrostatic pressure
  INIT_PRESSURE=0.0_CMISSRP !2.0_CMISSRP*C(2)-C(1)
  CALL cmfe_Field_ParametersToFieldParametersComponentCopy(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 1,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,Err)
  CALL cmfe_Field_ParametersToFieldParametersComponentCopy(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 2,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,2,Err)
  CALL cmfe_Field_ParametersToFieldParametersComponentCopy(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 3,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,3,Err)
  CALL cmfe_Field_ComponentValuesInitialise(DependentField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,4,INIT_PRESSURE,&
    & Err)
  !CALL cmfe_Field_ParametersToFieldParametersComponentCopy(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
  !  & 1,DependentField,CMFE_FIELD_V_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,Err)
  !CALL cmfe_Field_ParametersToFieldParametersComponentCopy(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
  !  & 2,DependentField,CMFE_FIELD_V_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,2,Err)
  !CALL cmfe_Field_ParametersToFieldParametersComponentCopy(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
  !  & 3,DependentField,CMFE_FIELD_V_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,3,Err)

!
! Problem and solver---------------------------------------------------------------------------------------------------------START
!

  ! Define the problem
  CALL cmfe_Problem_Initialise(Problem,Err)
  CALL cmfe_Problem_CreateStart(ProblemUserNumber,[CMFE_PROBLEM_ELASTICITY_CLASS,CMFE_PROBLEM_FINITE_ELASTICITY_TYPE, &
    & CMFE_PROBLEM_NO_SUBTYPE],Problem,Err)
  CALL cmfe_Problem_CreateFinish(Problem,Err)

  ! Create the problem control loop
  CALL cmfe_Problem_ControlLoopCreateStart(Problem,Err)
  !CALL cmfe_ControlLoop_In itialise(ControlLoop,Err)
  !CALL cmfe_Problem_ControlLoopGet(Problem,CMFE_CONTROL_LOOP_NODE,ControlLoop,Err)
  !CALL cmfe_ControlLoop_TypeSet(ControlLoop,CMFE_PROBLEM_CONTROL_LOAD_INCREMENT_LOOP_TYPE,Err)
  !CALL cmfe_ControlLoop_MaximumIterationsSet(ControlLoop,10,Err)
  !CALL cmfe_ControlLoop_LoadOutputSet(ControlLoop,1,Err)
  CALL cmfe_Problem_ControlLoopCreateFinish(Problem,Err)

  !Create the problem solvers
  CALL cmfe_Solver_Initialise(Solver,Err)
  CALL cmfe_Solver_Initialise(MechanicalSolver,Err)
  CALL cmfe_Problem_SolversCreateStart(Problem,Err)
  CALL cmfe_Problem_SolverGet(Problem,CMFE_CONTROL_LOOP_NODE,1,Solver,Err)
  CALL cmfe_Solver_OutputTypeSet(Solver,CMFE_SOLVER_PROGRESS_OUTPUT,Err)
  CALL cmfe_Solver_NewtonJacobianCalculationTypeSet(Solver,CMFE_SOLVER_NEWTON_JACOBIAN_EQUATIONS_CALCULATED,Err)
  CALL cmfe_Solver_NewtonLinearSolverGet(Solver,MechanicalSolver,Err)
  !CALL cmfe_Solver_LinearTypeSet(MechanicalSolver,CMFE_SOLVER_LINEAR_DIRECT_SOLVE_TYPE,Err)
  CALL cmfe_Solver_LinearTypeSet(MechanicalSolver,CMFE_SOLVER_NONLINEAR_NEWTON,Err)  
  CALL cmfe_Solver_LibraryTypeSet(MechanicalSolver,CMFE_SOLVER_MUMPS_LIBRARY,Err)
  CALL cmfe_Solver_NewtonRelativeToleranceSet(Solver,1.E-6_CMISSRP,Err)
  !CALL cmfe_Solver_NewtonAbsoluteToleranceSet(Solver,1.E-6_CMISSRP,Err)
  CALL cmfe_Solver_NewtonMaximumIterationsSet(Solver,300,Err)
  CALL cmfe_Problem_SolversCreateFinish(Problem,Err)

  !Create the problem solver equations
  CALL cmfe_Solver_Initialise(Solver,Err)
  CALL cmfe_SolverEquations_Initialise(SolverEquations,Err)
  CALL cmfe_Problem_SolverEquationsCreateStart(Problem,Err)
  CALL cmfe_Problem_SolverGet(Problem,CMFE_CONTROL_LOOP_NODE,1,Solver,Err)
  CALL cmfe_Solver_SolverEquationsGet(Solver,SolverEquations,Err)
  CALL cmfe_SolverEquations_EquationsSetAdd(SolverEquations,EquationsSet,EquationsSetIndex,Err)
  CALL cmfe_Problem_SolverEquationsCreateFinish(Problem,Err)

!
! Boundary conditions---------------------------------------------------------------------------------------------------------START
!

  !Prescribe boundary conditions (absolute nodal parameters)
  CALL cmfe_BoundaryConditions_Initialise(BoundaryConditions,Err)
  CALL cmfe_SolverEquations_BoundaryConditionsCreateStart(SolverEquations,BoundaryConditions,Err)

  !CALL cmfe_GeneratedMesh_SurfaceGet(GeneratedMesh,CMFE_GENERATED_MESH_REGULAR_BOTTOM_SURFACE,BottomSurfaceNodes,BottomNormalXi, &
    !& Err)
  !CALL cmfe_GeneratedMesh_SurfaceGet(GeneratedMesh,CMFE_GENERATED_MESH_REGULAR_LEFT_SURFACE,LeftSurfaceNodes,LeftNormalXi,Err)
  !CALL cmfe_GeneratedMesh_SurfaceGet(GeneratedMesh,CMFE_GENERATED_MESH_REGULAR_RIGHT_SURFACE,RightSurfaceNodes,RightNormalXi,Err)
  !CALL cmfe_GeneratedMesh_SurfaceGet(GeneratedMesh,CMFE_GENERATED_MESH_REGULAR_FRONT_SURFACE,FrontSurfaceNodes,FrontNormalXi,Err)
  !CALL cmfe_GeneratedMesh_SurfaceGet(GeneratedMesh,CMFE_GENERATED_MESH_REGULAR_TOP_SURFACE,TopSurfaceNodes,RightNormalXi,Err)
  !CALL cmfe_GeneratedMesh_SurfaceGet(GeneratedMesh,CMFE_GENERATED_MESH_REGULAR_BACK_SURFACE,BackSurfaceNodes,RightNormalXi,Err)

  !Set x=0 nodes to no x displacment in x. Set x=WIDTH nodes to 10% x displacement
  !DO node_idx=1,SIZE(LeftSurfaceNodes,1)
    !NodeNumber=LeftSurfaceNodes(node_idx)
    !CALL cmfe_Decomposition_NodeDomainGet(Decomposition,NodeNumber,1,NodeDomain,Err)
    !IF(NodeDomain==ComputationalNodeNumber) THEN
      !CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,NodeNumber,1, &
        !& CMFE_BOUNDARY_CONDITION_FIXED,0.0_CMISSRP,Err)
    !ENDIF
  !ENDDO


  !DO node_idx=1,SIZE(RightSurfaceNodes,1)
    !NodeNumber=RightSurfaceNodes(node_idx)
    !CALL cmfe_Decomposition_NodeDomainGet(Decomposition,NodeNumber,1,NodeDomain,Err)
    !IF(NodeDomain==ComputationalNodeNumber) THEN
      !CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,NodeNumber,1, &
        !& CMFE_BOUNDARY_CONDITION_FIXED,1.0_CMISSRP*WIDTH,Err)
    !ENDIF
  !ENDDO

  !Set y=0 nodes to no y displacement
  !DO node_idx=1,SIZE(FrontSurfaceNodes,1)
    !NodeNumber=FrontSurfaceNodes(node_idx)
    !CALL cmfe_Decomposition_NodeDomainGet(Decomposition,NodeNumber,1,NodeDomain,Err)
    !IF(NodeDomain==ComputationalNodeNumber) THEN
      !CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,NodeNumber,2, &
        !& CMFE_BOUNDARY_CONDITION_FIXED,0.0_CMISSRP,Err)
    !ENDIF
  !ENDDO

  !Set z=0 nodes to no z displacement
  !DO node_idx=1,SIZE(BottomSurfaceNodes,1)
    !NodeNumber=BottomSurfaceNodes(node_idx)
    !CALL cmfe_Decomposition_NodeDomainGet(Decomposition,NodeNumber,1,NodeDomain,Err)
    !IF(NodeDomain==ComputationalNodeNumber) THEN
      !CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,NodeNumber,3, &
        !& CMFE_BOUNDARY_CONDITION_FIXED,0.0_CMISSRP,Err)
    !ENDIF
  !ENDDO

  !Fix the top tendon nodes in all directions
  ! removed _INCREMENTED from BC type 
  VALUE = 0.0_CMISSRP
  DO node_idx=1,SIZE(TendonTopFixNodes,1)
    NodeNumber=TendonTopFixNodes(node_idx)
    CALL cmfe_Decomposition_NodeDomainGet(Decomposition,NodeNumber,1,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN

      CALL cmfe_Field_ParameterSetGetNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,NodeNumber,& 
        & 1,VALUE,Err)
      !VALUE=VALUE-15.8_CMISSRP
      CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,NodeNumber,1, &
        & CMFE_BOUNDARY_CONDITION_FIXED,VALUE,Err)

      CALL cmfe_Field_ParameterSetGetNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,NodeNumber,& 
        & 2,VALUE,Err)
      !VALUE=VALUE+5.2_CMISSRP
      CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,NodeNumber,2, &
        & CMFE_BOUNDARY_CONDITION_FIXED,VALUE,Err)

      CALL cmfe_Field_ParameterSetGetNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,NodeNumber,& 
        & 3,VALUE,Err)
      !VALUE=VALUE+31.8_CMISSRP
      CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,NodeNumber,3, &
        & CMFE_BOUNDARY_CONDITION_FIXED,VALUE,Err)
    ENDIF
  ENDDO

  !Fix the bottom tendon nodes in all directions
  DO node_idx=1,SIZE(TendonBottomFixNodes,1)
    NodeNumber=TendonBottomFixNodes(node_idx)
    CALL cmfe_Decomposition_NodeDomainGet(Decomposition,NodeNumber,1,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN

      CALL cmfe_Field_ParameterSetGetNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,NodeNumber,& 
        & 1,VALUE,Err)
      !VALUE=VALUE+10.0_CMISSRP
      CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,NodeNumber,1, &
        & CMFE_BOUNDARY_CONDITION_FIXED,VALUE,Err)

      CALL cmfe_Field_ParameterSetGetNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,NodeNumber,& 
        & 2,VALUE,Err)
      !VALUE=VALUE+1.0_CMISSRP
      CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,NodeNumber,2, &
        & CMFE_BOUNDARY_CONDITION_FIXED,VALUE,Err)

      CALL cmfe_Field_ParameterSetGetNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,NodeNumber,& 
        & 3,VALUE,Err)
      !VALUE=VALUE-17.7_CMISSRP
      CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,NodeNumber,3, &
        & CMFE_BOUNDARY_CONDITION_FIXED,VALUE,Err)
    ENDIF
  ENDDO


  CALL cmfe_SolverEquations_BoundaryConditionsCreateFinish(SolverEquations,Err)

!
! Solving the problem---------------------------------------------------------------------------------------------------------START
!

  CALL cmfe_Fields_Initialise(Fields,Err)
  CALL cmfe_Fields_Create(Region,Fields,Err)

  ! Compute the distances from the activation (sEMG) sources and each SIM-FIELD point-----------------------------------------START
  ! Where SIM-FIELD point = {NODE point, GAUSS point, ...} 

  
  ! 1) Define (read in) activation source(s): sEMG(j) = {EMGX_j, EMGY_j, EMGZ_j}, j = 1,...,EMG_TOT
  ! >>> currently only consider a single activation source >>>
    !activation_Source = (/0.5_CMISSRP,0.0_CMISSRP,0.5_CMISSRP/)

  ! read in activation sources
    OPEN(UNIT=12, FILE="activation_source_loc.txt", STATUS='OLD', ACTION='read')
      READ(12,*) activation_Source
    CLOSE(12)  
    ! note that the read in is transposed, so in the text file: EMG(1,1), EMGY(1,2), EMG(1,3) ; EMG(2,1), EMG(2,2), EMG(2,3)
    ! then the array will be a_S(1,1) = EMG(1,1), a_S(2,1) = EMG(1,2), a_S(3,1) = EMG(1,3) and so on 

  ! AT TIME t=t_n, t = 0,..., t_TOT
  ! >>> currently consider no time evolution (only based on inital configuration) >>>

  !   2) Read in SIM-FIELD points - here node points - Get nodal coordinates of all node points 
  !   NC_i(t_n) = {NCX_i(t_n), NCY_i(t_n), NCZ_i(t_n)}, i = 1,...,N_TOT 
      nodal_Coords = (/0.0_CMISSRP,0.0_CMISSRP,0.0_CMISSRP/)
      allocate(dist_Node_actSource(TotalNumberOfNodes, TotalNumberOfSources))
      allocate(alpha_spatial_wt(TotalNumberOfNodes, TotalNumberOfSources))
      DO node_idx=1,TotalNumberOfNodes
        CALL cmfe_Decomposition_NodeDomainGet(Decomposition,node_idx,1,NodeDomain,Err)
        IF(NodeDomain==ComputationalNodeNumber) THEN      
          DO component_idx=1,3   
            CALL cmfe_Field_ParameterSetGetNode(RegionUserNumber,FieldDependentUserNumber,CMFE_FIELD_U_VARIABLE_TYPE, &
              & CMFE_FIELD_VALUES_SET_TYPE, 1,1,node_idx,component_idx,nodal_Coords(component_idx),err)
          ENDDO
        
          DO source_idx=1,TotalNumberOfSources
  !         3) Compute the distance between each activation source and all nodal points
  !           d_ij(t_n) = sqrt { [EMGX_j-NCX_i(t_n)]² + [EMGY_j-NCY_i(t_n)]² + [EMGZ_j-NCZ_i(t_n)]² }
  !           d_ij is a 2D matrix with dimensions EMG_TOTxN_TOT for all time steps t_TOT  
            dist_Node_actSource(node_idx, source_idx) = sqrt( (nodal_Coords(1)-activation_Source(1,source_idx))**2 & 
              & + (nodal_Coords(2)-activation_Source(2,source_idx))**2 + (nodal_Coords(3)-activation_Source(3,source_idx))**2)
  !         4) Evaluate the spatial activation weighting function for each activation source, i.e.
  !           alpha_spatial_wt(NC_i(t_n)) = f(d_ij(t_n)) E [0,1], e.g.
  !           f(d_ij(t_n)) = exp( d_ij(t_n) * ( ln(gamma)/beta ) ),
  !             where beta is the maximum distance of influence and gamma is "some small value" 
            alpha_spatial_wt(node_idx, source_idx) = exp( dist_Node_actSource(node_idx, source_idx) * log(gamma)/beta )
  !           i) All the contributions need to be summed, i.e. the contribution of all activation sources to the current node
  !             alpha_spatial_wt_sum(NC_i(t_n)) = sum (f(d_ij(t_n))) for a given node over all points
  !           ii) The saturation activation is 1 
          ENDDO
        ENDIF
  !   5) Multiply the spatial activation weighting function with the temporal activation, i.e.
  !     alpha(NC_i(t_n)) = alpha_spatial_wt(NC_i(t_n)) * alpha(t_n)
      ENDDO 
  !   6) Now the spatial and temporal activation value can be passed into the muslce model for that given node
  !     F_active= f(alpha(NC_i(t_n)), ...)

  WRITE(*,*) "Completed computing spatial activations..."

  ! Spatial activation ---------------------------------------------------------------------------------------------------------END


  ! Read in activation alpha at time t; no check currently made that TIMESTEP = number of rows in text file
  OPEN (UNIT=3, FILE='activation_alpha.txt', STATUS='OLD', ACTION='read')
  DO i=1,TIMESTEPS
    READ(3,*) alpha(i)
  ENDDO
  CLOSE(3)
  
  ! Solve the mechanical problem fpr this time step - DUMMY SOLVE - required to initialise all variables 
  ! (ideally should not be required)
  CALL cmfe_Problem_Solve(Problem,Err)

  ! Mechanical BC increment (if concentric contraction is used)
  BCLOAD = 0.0_CMISSRP

  ! Loop over Time -----------------------------------------------------------------------------------------------------------START
  DO i=1,TIMESTEPS
    WRITE(*,*) "----------------------------------------------------------------------"
    WRITE(*,*) "Applying temporal activation, STEP: ", i
    ! Initialise alpha(i=t) 
    !CALL cmfe_Field_ComponentValuesInitialise(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,5, &
      !& 0.0_CMISSRP,Err)

    ! Set how alpha evolves over time (may choose to keep it constant, scale it, etc.)
    alpha_t = 0.0125_CMISSRP*alpha(i)
    WRITE(*,*) "Temporal Activation: ", alpha_t
    WRITE(*,*) "----------------------------------------------------------------------"
    ! Below are different options to update alpha at each time step ----------------------------------------------------------START
      ! (1) Entire muscle (working)
      ! (2) Element based (working - need to check consistancy with (1) & (3)): currently only valid option for spatial variation
      ! (3) Node based (working - need to check consistancy with (1) & (2))
      ! (4) Gauss-point based (not yet tested)

      ! Set spatial update case - WARNING: cmfe_Field_ComponentInterpolationSet for C(5) needs to be correspondingly set!  
      alpha_SU = 2

      SELECT CASE (alpha_SU)
        CASE(1) ! Entire muscle (working)
          CALL cmfe_Field_ComponentValuesInitialise(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,9, &
            & alpha_t,Err)

        CASE(2) ! Element based (working)
          ! loop over all elements
          DO elem_idx=1,23
            CALL cmfe_Field_ParameterSetUpdateElement(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,&
              & elem_idx,9,alpha_t,Err)
          END DO        
          DO elem_idx=28,NumberOfElements
            CALL cmfe_Field_ParameterSetUpdateElement(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,&
            & elem_idx,9, alpha_t,Err)
          END DO        
        CASE(3) ! Node based
          
          ! loop over all nodes 
          DO node_idx=1,TotalNumberOfNodes 
            CALL cmfe_Decomposition_NodeDomainGet(Decomposition,node_idx,1,NodeDomain,Err)
            IF(NodeDomain==ComputationalNodeNumber) THEN                
              alpha_spatial_wt_sum = 0.0_CMISSRP
              ! Apply spatial weighting according to actiation source
              DO source_idx=1,TotalNumberOfSources
                alpha_spatial_wt_sum =  alpha_spatial_wt_sum + alpha_spatial_wt(node_idx,source_idx)
              ENDDO
              ! Saturation of activation at alpha = 1
              IF (alpha_spatial_wt_sum .GE. 1.0_CMISSRP) THEN
                alpha_spatial_wt_sum = 1.0_CMISSRP
              ENDIF
              alpha_spatial_temporal = alpha_t * alpha_spatial_wt_sum
              write(*,*) alpha_spatial_temporal
              CALL cmfe_Field_ParameterSetUpdateNode(MaterialField,& 
                & CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1, 1, node_idx,9,alpha_spatial_temporal,Err)      
            ENDIF       
          ENDDO   

        CASE DEFAULT ! Entire muscle (working)
           CALL cmfe_Field_ComponentValuesInitialise(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,5, &
            & alpha_t,Err)     
      END SELECT! Spatial alpha update------------------------------------------------------------------------------------------END

    ! Solve the mechanical problem fpr this time step
    CALL cmfe_Problem_Solve(Problem,Err)
  
    ! update the mechanical boundary condition -------------------------------------------------------------------------------START
      !DO node_idx=1,SIZE(RightSurfaceNodes,1)
        !NodeNumber=RightSurfaceNodes(node_idx)
        !CALL cmfe_Decomposition_NodeDomainGet(Decomposition,NodeNumber,1,NodeDomain,Err)
        !IF(NodeDomain==ComputationalNodeNumber) THEN
          !CALL cmfe_Field_ParameterSetAddNode(DependentField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, & 
            !& 1,1,NodeNumber,1, BCLOAD,Err)
        !ENDIF
      !ENDDO ! mechanical BC update ---------------------------------------------------------------------------------------------END
 
    ! Output the exnode file at current time step 
    ! NOTE: if changing filename remember to change AXY!
    IF (i.LE.9) THEN
      WRITE(filename, "(A22,I1)") "ActiveStrain_BicepsExample_", i
      filename=trim(filename)
    ELSEIF (i.LE.99) THEN
      WRITE(filename, "(A22,I2)") "ActiveStrain_BicepsExample_", i
      filename=trim(filename)
    ELSE  
      WRITE(filename, "(A22,I3)") "ActiveStrain_BicepsExample_", i
      filename=trim(filename)
    ENDIF
    CALL cmfe_Fields_NodesExport(Fields,filename,"FORTRAN",Err)
  
  END DO
  ! Loop over Time -------------------------------------------------------------------------------------------------------------END  
  
  !Output solution
  CALL cmfe_Fields_NodesExport(Fields,"ActiveStrain_BicepsExample","FORTRAN",Err)
  CALL cmfe_Fields_ElementsExport(Fields,"ActiveStrain_BicepsExample","FORTRAN",Err)
  CALL cmfe_Fields_Finalise(Fields,Err)

!
! ---------------------------------------------------------------------------------------------------------------------------------
!

  CALL cmfe_Finalise(Err)

  WRITE(*,'(A)') "Program successfully completed."

  STOP

END PROGRAM LARGEUNIAXIALEXTENSIONEXAMPLE



