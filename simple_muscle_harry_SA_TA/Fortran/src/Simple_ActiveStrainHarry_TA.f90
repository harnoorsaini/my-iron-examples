! ---------------------------------------------------------------------------------------------------------------------------------
!
! SIMULATION OF A 3D HILL TYPE MUSCLE MODEL
!
! - Displacement controlled - isometric or concentric
! - Muscle activation - variation with time (read in from text file)
! - Muscle activation - variation over space
! - User (unstructed) mesh
!
! AUTHOR: HARNOOR SAINI
! MARCH 2017
!
! +++++++++++++++++++++++++++++++++++++++++++++++++++
! TO DO
! [ ] Check unit system - must be consistant: geometr & material parameters
! [x] Change material model to EQUATIONS_SET_TRANS_ISOTROPIC_ACTIVE_TRANSITION_SUBTYPE (try in simple model first)
! [ ] What about the pressure basis and incompressiblity...
!     --> Simple example does not converge when UserPressureBasis is off...
! [x] How are the fibre directions actually defined - it does not seem wrt global coordinates rather in the xi1 direction or 
!     something - leave the transforamtion to 0 for plausible fibre orientation (based on the specific node numbering)
! [ ] Where do the constants in the material models come from, e.g. fibre length bounds for active part @ line 4541 in finite_
!     elasticity_routines.f90?
! [x] Set up time loop as in TA example... (been there, done that - no advantage really)
! [ ] Why are the muslce forces flipping sings (flip!)
! [ ] What about the initial pressure value? - ask Andreas
! +++++++++++++++++++++++++++++++++++++++++++++++++++
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
  ! 
  ! UNITS
  ! mm, N, MPa, tonne
  ! source for Rivlin-skin: The uniaxial stress versus strain response of pig skin
  !                          and silicone rubber at low and high strain rates: C1: 0.3MPa, C2: 0MPa
  ! source for Rivlin-muscle: SPrenger PhD
  ! source for Markert-skin: see above paper - collagen Di
  ! source for Markert-muscle: see above paper - collagen VLi
  !REAL(CMISSRP), PARAMETER, DIMENSION(11) :: C= &
    !& [0.356_CMISSRP,0.386_CMISSRP,0.3411E-3_CMISSRP, &
    !&  44.0_CMISSRP,0.3_CMISSRP, 0.01_CMISSRP, 0.00646E-3_CMISSRP, &
    !&  30.0_CMISSRP, 0.0_CMISSRP, 1.0_CMISSRP, 0.5_CMISSRP] 

  ! "original" parameters
  REAL(CMISSRP), PARAMETER, DIMENSION(11) :: C= &
    & [3.56E-2_CMISSRP,3.86E-2_CMISSRP,0.3E-8_CMISSRP, &
    &  34.0_CMISSRP,3.56E-2_CMISSRP, 3.86E-2_CMISSRP, 0.3E-8_CMISSRP, &
    &  34.0_CMISSRP, 0.0_CMISSRP, 1.0_CMISSRP, 0.5_CMISSRP]   

  ! Test program parameters
  REAL(CMISSRP), PARAMETER :: MPc = 1.0_CMISSRP ! Material Parameter constant multiplied; ONLY PARAMETERS 1-8
  REAL(CMISSDP), PARAMETER :: PI=4.0_CMISSRP*DATAN(1.0_CMISSRP)
  REAL(CMISSRP), PARAMETER :: PERIOD=1.0_CMISSRP
  REAL(CMISSRP), PARAMETER :: TIME_STOP =10000.0_CMISSRP
  REAL(CMISSRP), PARAMETER :: alpha_inc = 1.0E-5_CMISSRP
  

  !INTEGER(CMISSIntg), PARAMETER :: InterpolationType=CMFE_BASIS_LINEAR_LAGRANGE_INTERPOLATION
  INTEGER(CMISSIntg), PARAMETER :: InterpolationType=CMFE_BASIS_QUADRATIC_LAGRANGE_INTERPOLATION
  INTEGER(CMISSIntg), PARAMETER :: PressureInterpolationType=CMFE_BASIS_LINEAR_LAGRANGE_INTERPOLATION
  LOGICAL, PARAMETER :: UsePressureBasis=.FALSE. ! equivalent to TA example...
  INTEGER(CMISSIntg), PARAMETER :: NumberOfGaussXi=3

  INTEGER(CMISSIntg), PARAMETER :: CoordinateSystemUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: RegionUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: BasisUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: PressureBasisUserNumber=3
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
  
  ! Program types
  ! Program variables
  INTEGER(CMISSIntg) :: NumberGlobalXElements,NumberGlobalYElements,NumberGlobalZElements
  INTEGER(CMISSIntg) :: EquationsSetIndex, TotalNumberOfNodes
  INTEGER(CMISSIntg) :: NumberOfComputationalNodes,NumberOfDomains,ComputationalNodeNumber
  INTEGER(CMISSIntg) :: NodeNumber,NodeDomain,node_idx,elem_idx, gauss_idx, component_idx, my_node_idx
  INTEGER(CMISSIntg),ALLOCATABLE :: BottomSurfaceNodes(:)
  INTEGER(CMISSIntg),ALLOCATABLE :: LeftSurfaceNodes(:)
  INTEGER(CMISSIntg),ALLOCATABLE :: RightSurfaceNodes(:)
  INTEGER(CMISSIntg),ALLOCATABLE :: FrontSurfaceNodes(:)
  INTEGER(CMISSIntg),ALLOCATABLE :: TopSurfaceNodes(:)
  INTEGER(CMISSIntg),ALLOCATABLE :: BackSurfaceNodes(:)
  INTEGER(CMISSIntg) :: BottomNormalXi,LeftNormalXi,RightNormalXi,FrontNormalXi
  INTEGER(CMISSIntg) :: FibreFieldNumberOfComponents
  !INTEGER(CMISSIntg), PARAMETER :: NUMBER_OF_COMPONENTS = 3 !nearly incompressible
  INTEGER(CMISSIntg), PARAMETER :: NUMBER_OF_COMPONENTS = 4 !fully incompressible

  ! CMISS variables
  TYPE(cmfe_ControlLoopType) :: ControlLoopMain
  TYPE(cmfe_BasisType) :: QuadraticBasis,LinearBasis, PressureBasis
  TYPE(cmfe_BoundaryConditionsType) :: BoundaryConditions
  TYPE(cmfe_CoordinateSystemType) :: CoordinateSystem, WorldCoordinateSystem
  TYPE(cmfe_MeshType) :: Mesh
  TYPE(cmfe_DecompositionType) :: Decomposition
  TYPE(cmfe_EquationsType) :: Equations
  TYPE(cmfe_EquationsSetType) :: EquationsSet
  TYPE(cmfe_FieldType) :: GeometricField,MaterialField,DependentField,EquationsSetField,FibreField
  TYPE(cmfe_FieldsType) :: Fields
  TYPE(cmfe_GeneratedMeshType) :: GeneratedMesh
  TYPE(cmfe_ProblemType) :: Problem
  TYPE(cmfe_RegionType) :: Region,WorldRegion
  TYPE(cmfe_SolverType) :: Solver,LinearSolver
  TYPE(cmfe_SolverEquationsType) :: SolverEquations
  TYPE(cmfe_ControlLoopType) :: ControlLoop
  TYPE(cmfe_NodesType) :: Nodes
  TYPE(cmfe_MeshElementsType) :: QuadraticElements
  TYPE(cmfe_MeshElementsType) :: LinearElements
  TYPE(cmfe_MeshElementsType) :: ElementsM

  INTEGER(CMISSIntg) :: i
  CHARACTER(LEN=256) :: filename
  
#ifdef WIN32
  !Quickwin type
  LOGICAL :: QUICKWIN_STATUS=.FALSE.
  TYPE(WINDOWCONFIG) :: QUICKWIN_WINDOW_CONFIG
#endif

  ! Generic CMISS variables
  INTEGER(CMISSIntg) :: Err

  ! Variables, Parameters, ...for specific simulation
  !INTEGER(CMISSIntg), PARAMETER :: TIMESTEPS=11 ! Number of Timesteps
  INTEGER(CMISSIntg) :: alpha_SU, j, k, stat ! spatial update of alpha
  !REAL(CMISSRP), DIMENSION(TIMESTEPS) :: alpha
  REAL(CMISSRP) :: alpha_t, BCLOAD, VALUE, time
  REAL(CMISSRP) :: FibreFieldAngle(3)
!
! User Geometric--------------------------------------------------------------------------------------------------------------START
!

! the nodes in each element for quadratic basis functions
  INTEGER(CMISSIntg), PARAMETER :: NumberOfElementsFE=39
  INTEGER(CMISSIntg), PARAMETER :: NumberOfNodesFE=543
  INTEGER(CMISSIntg), DIMENSION(NumberOfElementsFE,27) :: AllQuadraticElements
  character(len=8), parameter :: element_str = "Element:"
  character(len=6), parameter :: nodes_str = "Nodes:"
  character(len=4), parameter :: node_str = "Node"
  character(len=5), parameter :: nodep_str = "Node:"
  character(len=132) :: blub,blub2,blub3,blub4,blub5,blub6

! For boundary conditions
  INTEGER(CMISSIntg), PARAMETER, DIMENSION(21) :: TOP_NODES = &
   & [66,67,127,206,213,217,221,254,255,262,263,290,291,297,299,340,341,344,345,347,363]
  INTEGER(CMISSIntg), PARAMETER, DIMENSION(21) :: BOTTOM_NODES = &
   & [57,58,107,229,230,231,232,232,234,281,282,283,284,326,327,328,329,330,331,360,361]
  INTEGER(CMISSIntg), PARAMETER, DIMENSION(52) :: SIDE_NODES = &
   & [484,485,500,220,222,261,486,487,501,214,251,265,488,489,502,216,218,267,490,491,503, &
   &  209,240,269,492,493,504,212,219,276,494,495,505,224,228,278,496,497,506,225,227,280, &
   &  538,540,542,233,498,499,507,539,541,543]

INTEGER(CMISSIntg), DIMENSION(NumberOfElementsFE,27), PARAMETER :: AllElementNodes=transpose(reshape( &
   & [6,21,103,31,33,155,8,22,104,46,47,182,48,49,184,70,71,230,13,23,125,32,34,156,15,24,127, &
   &  13,23,125,32,34,156,15,24,127,54,55,202,56,57,204,74,75,240,9,25,121,35,37,169,11,26,123, &
   &  9,25,121,35,37,169,11,26,123,50,51,198,52,53,200,72,73,238,17,27,129,36,38,170,19,28,131, &
   &  17,27,129,36,38,170,19,28,131,58,59,206,60,61,208,76,77,242,120,194,101,166,196,165,117,236,102, &
   &  103,137,105,155,157,173,104,138,106,182,183,214,184,185,215,230,231,246,125,139,126,156,158,174,127,140,128, &
   &  109,141,113,159,161,175,110,142,114,186,187,216,188,189,217,232,233,247,133,143,134,160,162,176,135,144,136, &
   &  111,145,115,163,164,177,112,146,116,190,191,218,192,193,219,234,235,248,103,137,105,155,157,173,104,138,106, &
   &  101,147,107,165,167,178,102,148,108,194,195,220,196,197,221,236,237,249,120,149,119,166,168,179,117,150,118, &
   &  121,151,122,169,171,180,123,152,124,198,199,222,200,201,223,238,239,250,129,153,130,170,172,181,131,154,132, &
   &  125,139,126,156,158,174,127,140,128,202,203,224,204,205,225,240,241,251,121,151,122,169,171,180,123,152,124, &
   &  129,153,130,170,172,181,131,154,132,206,207,226,208,209,227,242,243,252,101,147,107,165,167,178,102,148,108, &
   &  133,143,134,160,162,176,135,144,136,210,211,228,212,213,229,244,245,253,111,145,115,163,164,177,112,146,116, &
   &  1017,1075,111,1110,1112,163,1068,1076,112,1190,1191,190,1192,1193,192,1238, &
   &  1239,234,1020,1073,103,1109,1111,155,1070,1074,104, &
   &  1051,1077,1001,1113,1115,1165,1033,1078,1035,1194,1195,1288,1196,1197,1289, &
   &  1334,1335,1350,1030,1079,1026,1114,1116,1166,1034,1080,1016, &
   &  1065,1081,1015,1117,1118,1167,1066,1082,1027,1198,1199,1290,1200,1201,1291, &
   &  1336,1337,1351,1051,1077,1001,1113,1115,1165,1033,1078,1035, &
   &  1059,1083,1009,1119,1120,1168,1031,1084,1029,1202,1203,1292,1204,1205,1293, &
   &  1338,1339,1352,1065,1081,1015,1117,1118,1167,1066,1082,1027, &
   &  1062,1085,1012,1121,1122,1169,1019,1086,1055,1206,1207,1294,1208,1209,1295, &
   &  1340,1341,1353,1059,1083,1009,1119,1120,1168,1031,1084,1029, &
   &  1058,1088,11,1124,1126,1171,1059,1083,1009,1210,1211,72, &
   &  1212,1213,1297,1202,1203,1292,1064,1087,19,1123,1125,1170,1065,1081,1015, &
   &  1020,1073,103,1109,1111,155,1070,1074,104,1214,1215,21,1216,1217,33, &
   &  1230,1231,22,1054,1089,6,1127,1128,31,1056,1090,8, &
   &  1061,1091,15,1129,1130,1173,1062,1085,1012,1218,1219,74, &
   &  1220,1221,1301,1206,1207,1294,1058,1088,11,1124,1126,1171,1059,1083,1009, &
   &  1056,1090,8,1131,1132,1174,1018,1092,1028,1222,1223,70, &
   &  1224,1225,1303,1226,1227,1304,1061,1091,15,1129,1130,1173,1062,1085,1012, &
   &  1018,1092,1028,1133,1134,1175,1032,1093,1025,1226,1227,1304, &
   &  1228,1229,1305,1342,1343,1354,1062,1085,1012,1121,1122,1169,1019,1086,1055, &
   &  1070,1074,104,1135,1136,1176,1053,1094,1024,1230,1231,22, &
   &  1232,1233,1307,1234,1235,1308,1056,1090,8,1131,1132,1174,1018,1092,1028, &
   &  1053,1094,1024,1137,1138,1177,1037,1095,1041,1234,1235,1308, &
   &  1236,1237,1309,1344,1345,1355,1018,1092,1028,1133,1134,1175,1032,1093,1025, &
   &  1068,1076,112,1139,1140,1178,1036,1096,1039,1238,1239,234, &
   &  1240,1241,1311,1242,1243,1312,1070,1074,104,1135,1136,1176,1053,1094,1024, &
   &  1036,1096,1039,1141,1142,1179,1040,1097,1038,1242,1243,1312, &
   &  1244,1245,1313,1346,1347,1356,1053,1094,1024,1137,1138,1177,1037,1095,1041, &
   &  1246,1247,133,1248,1249,160,1250,1251,135,1358,1360,210, &
   &  1364,1366,212,1370,1372,244,1017,1075,111,1110,1112,163,1068,1076,112, &
   &  1250,1251,135,1252,1253,1317,1254,1255,1318,1370,1372,244, &
   &  1376,1378,1380,1382,1384,1386,1068,1076,112,1139,1140,1178,1036,1096,1039, &
   &  1254,1255,1318,1256,1257,1319,1348,1349,1357,1382,1384,1386, &
   &  1388,1390,1392,1394,1396,1398,1036,1096,1039,1141,1142,1179,1040,1097,1038, &
   &  1054,1089,6,1127,1128,31,1056,1090,8,1258,1259,46,1260,1261,48, &
   &  1222,1223,70,1060,1102,13,1149,1150,32,1061,1091,15, &
   &  1060,1102,13,1149,1150,32,1061,1091,15,1262,1263,54,1264,1265,56, &
   &  1218,1219,74,1057,1103,9,1151,1152,35,1058,1088,11, &
   &  1057,1103,9,1151,1152,35,1058,1088,11,1266,1267,50,1268,1269,52, &
   &  1210,1211,72,1063,1104,17,1153,1154,36,1064,1087,19, &
   &  1063,1104,17,1153,1154,36,1064,1087,19,1270,1271,58,1272,1273,60, &
   &  1278,1279,76,1050,1105,120,1155,1156,166,1052,1106,117, &
   &  1050,1105,120,1155,1156,166,1052,1106,117,1274,1275,149,1276,1277,168, &
   &  1282,1283,150,1072,1107,119,1157,1158,179,1071,1108,118, &
   &  1064,1087,19,1123,1125,1170,1065,1081,1015,1278,1279,76,1280,1281,1331, &
   &  1198,1199,1290,1052,1106,117,1159,1160,1188,1051,1077,1001, &
   &  1052,1106,117,1159,1160,1188,1051,1077,1001,1282,1283,150,1284,1285,1333, &
   &  1194,1195,1288,1071,1108,118,1161,1162,1189,1030,1079,1026, &
   &  1043,1098,109,1143,1144,159,1049,1099,110,1359,1361,186, &
   &  1365,1367,188,1371,1373,232,1246,1247,133,1248,1249,160,1250,1251,135, &
   &  1049,1099,110,1145,1146,1181,1044,1100,1046,1371,1373,232, &
   &  1377,1379,1381,1383,1385,1387,1250,1251,135,1252,1253,1317,1254,1255,1318, &
   &  1044,1100,1046,1147,1148,1182,1047,1101,1045,1383,1385,1387, &
   &  1389,1391,1393,1395,1397,1399,1254,1255,1318,1256,1257,1319,1348,1349,1357], &
   & [27,NumberOfElementsFE]))

 INTEGER(CMISSIntg), DIMENSION(NumberOfNodesFE), PARAMETER :: AllNodes= &
   & [6,8,9,11,13,15,17,19,21,22,23,24,25,26,27,28,31,32,33,34,35,36,37,38,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,70,71, &
   & 72,73,74,75,76,77,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,125,126, &
   & 127,128,129,130,131,132,133,134,135,136,137,138,139,140,141,142,143,144,145,146,147,148,149,150,151,152,153,154,155,156,157, &
   & 158,159,160,161,162,163,164,165,166,167,168,169,170,171,172,173,174,175,176,177,178,179,180,181,182,183,184,185,186,187,188, &
   & 189,190,191,192,193,194,195,196,197,198,199,200,201,202,203,204,205,206,207,208,209,210,211,212,213,214,215,216,217,218,219, &
   & 220,221,222,223,224,225,226,227,228,229,230,231,232,233,234,235,236,237,238,239,240,241,242,243,244,245,246,247,248,249,250, &
   & 251,252,253,1001,1009,1012,1015,1016,1017,1018,1019,1020,1024,1025,1026,1027,1028,1029,1030,1031,1032,1033,1034,1035,1036, &
   & 1037,1038,1039,1040,1041,1043,1044,1045,1046,1047,1049,1050,1051,1052,1053,1054,1055,1056,1057,1058,1059,1060,1061,1062, &
   & 1063,1064,1065,1066,1068,1070,1071,1072,1073,1074,1075,1076,1077,1078,1079,1080,1081,1082,1083,1084,1085,1086,1087,1088, &
   & 1089,1090,1091,1092,1093,1094,1095,1096,1097,1098,1099,1100,1101,1102,1103,1104,1105,1106,1107,1108,1109,1110,1111,1112, &
   & 1113,1114,1115,1116,1117,1118,1119,1120,1121,1122,1123,1124,1125,1126,1127,1128,1129,1130,1131,1132,1133,1134,1135,1136,1137,&
   & 1138,1139,1140,1141,1142,1143,1144,1145,1146,1147,1148,1149,1150,1151,1152,1153,1154,1155,1156,1157,1158,1159,1160,1161,1162,&
   & 1165,1166,1167,1168,1169,1170,1171,1173,1174,1175,1176,1177,1178,1179,1181,1182,1188,1189,1190,1191,1192,1193,1194,1195,1196,&
   & 1197,1198,1199,1200,1201,1202,1203,1204,1205,1206,1207,1208,1209,1210,1211,1212,1213,1214,1215,1216,1217,1218,1219,1220,1221,&
   & 1222,1223,1224,1225,1226,1227,1228,1229,1230,1231,1232,1233,1234,1235,1236,1237,1238,1239,1240,1241,1242,1243,1244,1245,1246,&
   & 1247,1248,1249,1250,1251,1252,1253,1254,1255,1256,1257,1258,1259,1260,1261,1262,1263,1264,1265,1266,1267,1268,1269,1270,1271,&
   & 1272,1273,1274,1275,1276,1277,1278,1279,1280,1281,1282,1283,1284,1285,1288,1289,1290,1291,1292,1293,1294,1295,1297,1301,1303,&
   & 1304,1305,1307,1308,1309,1311,1312,1313,1317,1318,1319,1331,1333,1334,1335,1336,1337,1338,1339,1340,1341,1342,1343,1344,1345,&
   & 1346,1347,1348,1349,1350,1351,1352,1353,1354,1355,1356,1357,1358,1359,1360,1361,1364,1365,1366,1367,1370,1371,1372,1373,1376,&
   & 1377,1378,1379,1380,1381,1382,1383,1384,1385,1386,1387,1388,1389,1390,1391,1392,1393,1394,1395,1396,1397,1398,1399]

  !the entries in MyElemen1,... to be used in linearElement1,...
  INTEGER(CMISSIntg), DIMENSION(8), PARAMETER :: ENTRIES = [1,3,7,9,19,21,25,27]
  ! this should be 1 !!!!!
  real(CMISSRP), parameter :: scalefactor = 1.0_CMISSRP
  INTEGER(CMISSIntg), PARAMETER :: NumberOfSpatialCoordinates=3  
  INTEGER(CMISSIntg), PARAMETER :: NumberOfMeshComponentsFE=2
  INTEGER(CMISSIntg), PARAMETER :: QuadraticMeshComponentNumber=1
  INTEGER(CMISSIntg), PARAMETER :: LinearMeshComponentNumber=2
  INTEGER(CMISSIntg), PARAMETER :: MonodomainMeshComponentNumber=1  
  INTEGER(CMISSIntg), PARAMETER :: QuadraticBasisUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: LinearBasisUserNumber=2
  INTEGER(CMISSIntg), PARAMETER :: NumberOfXiCoordinates=NumberOfSpatialCoordinates
  INTEGER(CMISSIntg), PARAMETER :: NumberOfGaussPoints=3
  INTEGER(CMISSIntg), PARAMETER :: NumberOfMeshDimensionsFE=3
  real(CMISSRP) :: posX,posY,posZ
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

  !NumberGlobalXElements=2
  !NumberGlobalYElements=2
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

  ! Define geometric basis
  !Create the bases
  !Define basis functions - tri-Quadratic Lagrange 
  CALL cmfe_Basis_Initialise(QuadraticBasis,Err)
  CALL cmfe_Basis_CreateStart(QuadraticBasisUserNumber,QuadraticBasis,Err)
  CALL cmfe_Basis_TypeSet(QuadraticBasis,CMFE_BASIS_LAGRANGE_HERMITE_TP_TYPE,Err)
  CALL cmfe_Basis_NumberOfXiSet(QuadraticBasis,NumberOfXiCoordinates,Err)
  CALL cmfe_Basis_InterpolationXiSet(QuadraticBasis,[CMFE_BASIS_QUADRATIC_LAGRANGE_INTERPOLATION, &
   & CMFE_BASIS_QUADRATIC_LAGRANGE_INTERPOLATION,CMFE_BASIS_QUADRATIC_LAGRANGE_INTERPOLATION],Err)
  CALL cmfe_Basis_QuadratureNumberOfGaussXiSet(QuadraticBasis, &
   & [NumberOfGaussPoints,NumberOfGaussPoints,NumberOfGaussPoints],Err)
  CALL cmfe_Basis_CreateFinish(QuadraticBasis,Err)

  !Define basis functions - tri-Linear Lagrange
  CALL cmfe_Basis_Initialise(LinearBasis,Err)
  CALL cmfe_Basis_CreateStart(LinearBasisUserNumber,LinearBasis,Err)
  CALL cmfe_Basis_TypeSet(LinearBasis,CMFE_BASIS_LAGRANGE_HERMITE_TP_TYPE,Err)
  CALL cmfe_Basis_NumberOfXiSet(LinearBasis,NumberOfXiCoordinates,Err)
  CALL cmfe_Basis_InterpolationXiSet(LinearBasis,[CMFE_BASIS_LINEAR_LAGRANGE_INTERPOLATION, &
   & CMFE_BASIS_LINEAR_LAGRANGE_INTERPOLATION,CMFE_BASIS_LINEAR_LAGRANGE_INTERPOLATION],Err)
  CALL cmfe_Basis_QuadratureNumberOfGaussXiSet(LinearBasis, &
   & [NumberOfGaussPoints,NumberOfGaussPoints,NumberOfGaussPoints],Err)
  CALL cmfe_Basis_CreateFinish(LinearBasis,Err)

  ! Define pressure basis
  IF(UsePressureBasis) THEN
    CALL cmfe_Basis_Initialise(PressureBasis,Err)
    CALL cmfe_Basis_CreateStart(PressureBasisUserNumber,PressureBasis,Err)
    SELECT CASE(PressureInterpolationType)
    CASE(1,2,3,4)
      CALL cmfe_Basis_TypeSet(PressureBasis,CMFE_BASIS_LAGRANGE_HERMITE_TP_TYPE,Err)
    CASE(7,8,9)
      CALL cmfe_Basis_TypeSet(PressureBasis,CMFE_BASIS_SIMPLEX_TYPE,Err)
    END SELECT
    IF(NumberGlobalZElements==0) THEN
      CALL cmfe_Basis_NumberOfXiSet(PressureBasis,2,Err)
      CALL cmfe_Basis_InterpolationXiSet(PressureBasis,[PressureInterpolationType,PressureInterpolationType],Err)
      IF(NumberOfGaussXi>0) THEN
        CALL cmfe_Basis_QuadratureNumberOfGaussXiSet(PressureBasis,[NumberOfGaussXi,NumberOfGaussXi],Err)
      ENDIF
    ELSE
      CALL cmfe_Basis_NumberOfXiSet(PressureBasis,3,Err)
      CALL cmfe_Basis_InterpolationXiSet(PressureBasis, &
        & [PressureInterpolationType,PressureInterpolationType,PressureInterpolationType],Err)
      IF(NumberOfGaussXi>0) THEN
        CALL cmfe_Basis_QuadratureNumberOfGaussXiSet(PressureBasis,[NumberOfGaussXi,NumberOfGaussXi,NumberOfGaussXi],Err)
      ENDIF
    ENDIF
    CALL cmfe_Basis_CreateFinish(PressureBasis,Err)
  ENDIF

!
! User Mesh+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++START
!
!Create a mesh with three-dimensional elements
  CALL cmfe_Mesh_Initialise(Mesh,Err)
  CALL cmfe_Mesh_CreateStart(MeshUserNumber,Region,NumberOfMeshDimensionsFE,Mesh,Err)
  CALL cmfe_Mesh_NumberOfComponentsSet(Mesh,NumberOfMeshComponentsFE,Err) 
  CALL cmfe_Mesh_NumberOfElementsSet(Mesh,NumberOfElementsFE,Err)  
  !Define nodes for the mesh
  CALL cmfe_Nodes_Initialise(Nodes,Err)
  CALL cmfe_Nodes_CreateStart(Region,NumberOfNodesFE,Nodes,Err)
  CALL cmfe_Nodes_CreateFinish(Nodes,Err)


  ! README: 
  !
  ! Arbitrary user numbers in OpenCMISS does not really work. It is easier to leave the user numbers 
  ! identical to the global numbers (default).
  ! A mapping from the original CMISS numbering to the OpenCMISS consecutive numbering is therefore required:
  !
  ! AllNodes contains all node numbers that occur in the original CMISS mesh. These numbers are not consecutive. 
  !     This is a problem in OpenCMISS.
  ! Element1,... contain the node numbers corresponding to each Finite Element in the original CMISS numbering.
  ! 
  ! quadraticElement1,... contain the node numbers corresponding to each Finite Element in the new OpenCMISS numbering (1-81).
  DO i=1,27
    DO k=1,NumberOfNodesFE
      DO j=1,NumberOfElementsFE
        IF(AllElementNodes(j,i)==AllNodes(k)) AllQuadraticElements(j,i)=k
      ENDDO
    ENDDO
  ENDDO

  CALL cmfe_MeshElements_Initialise(QuadraticElements,Err)
  CALL cmfe_MeshElements_CreateStart(Mesh,QuadraticMeshComponentNumber,QuadraticBasis,QuadraticElements,Err)
  DO j=1,NumberOfElementsFE
    CALL cmfe_MeshElements_NodesSet(QuadraticElements,j,AllQuadraticElements(j,:),Err)
  ENDDO
  CALL cmfe_MeshElements_CreateFinish(QuadraticElements,Err)

  !for the linear elements, we need only specific entries from the elements above
  CALL cmfe_MeshElements_Initialise(LinearElements,Err)
  CALL cmfe_MeshElements_CreateStart(Mesh,LinearMeshComponentNumber,LinearBasis,LinearElements,Err)
  DO j=1,NumberOfElementsFE
    CALL cmfe_MeshElements_NodesSet(LinearElements,j,AllQuadraticElements(j,ENTRIES),Err)
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

!
! Decomposition---------------------------------------------------------------------------------------------------------------START
!

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
  ! From TA Example                                >>>> START  
  CALL cmfe_Field_TypeSet(GeometricField,CMFE_FIELD_GEOMETRIC_TYPE,Err)
  CALL cmfe_Field_MeshDecompositionSet(GeometricField,Decomposition,Err)  
  CALL cmfe_Field_NumberOfVariablesSet(GeometricField,1,Err)
  CALL cmfe_Field_NumberOfComponentsSet(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,3,Err) 
  CALL cmfe_Field_ComponentMeshComponentSet(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,1,QuadraticMeshComponentNumber,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,2,QuadraticMeshComponentNumber,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,3,QuadraticMeshComponentNumber,Err) 
  !CALL cmfe_Field_ScalingTypeSet(GeometricField,CMFE_FIELD_ARITHMETIC_MEAN_SCALING,Err)  
  !                                                >>>> END
  CALL cmfe_Field_VariableLabelSet(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,"Geometry",Err)
  CALL cmfe_Field_CreateFinish(GeometricField,Err)


!--------------------------------------------------------------------------------------------------------------------------------
  !Update the geometric field parameters for the finite elasticity geometric field 
  open(unit=2,file="input/FinalModell.ipnode",iostat=stat)
  do
    read(2,*,iostat=stat) blub,blub2,blub3,blub4,node_idx
    if(stat<0) exit !end of file
    !search for the lines that start with "Node" and store the node number
    if(blub==node_str) then
      !read the Xj(1), Xj(2), and Xj(3) position of the node
      read(2,*,iostat=stat) blub,blub2,blub3,blub4,blub5,blub6,posX
      if(stat/=0) write(*,*) "ERROR!!!"
      read(2,*,iostat=stat) blub,blub2,blub3,blub4,blub5,blub6,posY
      if(stat/=0) write(*,*) "ERROR!!!"
      read(2,*,iostat=stat) blub,blub2,blub3,blub4,blub5,blub6,posZ
      if(stat/=0) write(*,*) "ERROR!!!"
      !find the index of this node in AllNodes
      do i=1,size(AllNodes)
        if(node_idx==AllNodes(i)) then
          my_node_idx=i
          exit
        endif
      enddo
      !print the node number and its position
      !write(*,*) my_node_idx,posX,posY,posZ
      CALL cmfe_Decomposition_NodeDomainGet(Decomposition,my_node_idx,1,NodeDomain,Err)
      IF(NodeDomain==ComputationalNodeNumber) THEN
        CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1, &
         & my_node_idx,1,posX/scalefactor,Err)
        CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1, &
         & my_node_idx,2,posY/scalefactor,Err)
        CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1, &
         & my_node_idx,3,posZ/scalefactor,Err)
      ENDIF
    endif
  enddo
  close(unit=2)
  write(*,*) "Finished reading file: input/FinalModell.ipnode"

  CALL cmfe_Field_ParameterSetUpdateStart(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,Err)
  CALL cmfe_Field_ParameterSetUpdateFinish(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,Err)

  ! Update the geometric field parameters
  !CALL cmfe_GeneratedMesh_GeometricParametersCalculate(GeneratedMesh,GeometricField,Err)

  ! Create a fibre field and attach it to the geometric field
  CALL cmfe_Field_Initialise(FibreField,Err)
  CALL cmfe_Field_CreateStart(FieldFibreUserNumber,Region,FibreField,Err)
  CALL cmfe_Field_MeshDecompositionSet(FibreField,Decomposition,Err)  
  ! From TA Example                                >>>> START  
  CALL cmfe_Field_TypeSet(FibreField,CMFE_FIELD_GEOMETRIC_TYPE,Err)
  CALL cmfe_Field_NumberOfVariablesSet(FibreField,1,Err)
  CALL cmfe_Field_NumberOfComponentsSet(FibreField,CMFE_FIELD_U_VARIABLE_TYPE,3,Err) 
  !                                                >>>> END
  CALL cmfe_Field_TypeSet(FibreField,CMFE_FIELD_FIBRE_TYPE,Err)
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
  CALL cmfe_Field_ComponentInterpolationSet(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,1,CMFE_FIELD_CONSTANT_INTERPOLATION,Err)
  CALL cmfe_Field_ComponentInterpolationSet(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,2,CMFE_FIELD_CONSTANT_INTERPOLATION,Err)
  CALL cmfe_Field_ComponentInterpolationSet(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,3,CMFE_FIELD_CONSTANT_INTERPOLATION,Err)
  CALL cmfe_Field_ComponentInterpolationSet(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,4,CMFE_FIELD_CONSTANT_INTERPOLATION,Err)
  CALL cmfe_Field_ComponentInterpolationSet(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,5,CMFE_FIELD_CONSTANT_INTERPOLATION,Err)
  CALL cmfe_Field_ComponentInterpolationSet(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,6,CMFE_FIELD_CONSTANT_INTERPOLATION,Err)
  CALL cmfe_Field_ComponentInterpolationSet(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,7,CMFE_FIELD_CONSTANT_INTERPOLATION,Err)
  CALL cmfe_Field_ComponentInterpolationSet(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,8,CMFE_FIELD_CONSTANT_INTERPOLATION,Err)

  ! Spatial variation of muscle activation
    ! Required since value of activation is different at each point (see line 3303 in opencmiss_iron.f90)
    
    ! Constant interpolation (default)
    CALL cmfe_Field_ComponentInterpolationSet(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,9,CMFE_FIELD_CONSTANT_INTERPOLATION,Err)

    ! Element based interpolation
    !CALL cmfe_Field_ComponentInterpolationSet(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,9,CMFE_FIELD_ELEMENT_BASED_INTERPOLATION,Err)
    
    ! Node based interpolation
    !CALL cmfe_Field_ComponentInterpolationSet(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,5,CMFE_FIELD_NODE_BASED_INTERPOLATION,Err)
    
    ! Gauss-point based interpolation
    !CALL cmfe_Field_ComponentInterpolationSet(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,5, & 
    !& CMFE_FIELD_GAUSS_POINT_BASED_INTERPOLATION,Err)
  
  ! Spatial variation of passive/active tissue
    ! Constant
    ! CALL cmfe_Field_ComponentInterpolationSet(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,10,CMFE_FIELD_CONSTANT_INTERPOLATION,Err)

    ! Element based interpolation
    CALL cmfe_Field_ComponentInterpolationSet(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,10, &
      & CMFE_FIELD_ELEMENT_BASED_INTERPOLATION,Err)

  CALL cmfe_Field_ComponentInterpolationSet(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,11,CMFE_FIELD_CONSTANT_INTERPOLATION,Err) 
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
  CALL cmfe_Field_NumberOfComponentsSet(DependentField,CMFE_FIELD_U_VARIABLE_TYPE,NUMBER_OF_COMPONENTS,Err)
  CALL cmfe_Field_NumberOfComponentsSet(DependentField,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,NUMBER_OF_COMPONENTS,Err)  
  ! From TA Example                                >>>> START
  CALL cmfe_Field_ComponentMeshComponentSet(DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,QuadraticMeshComponentNumber,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(DependentField,CMFE_FIELD_U_VARIABLE_TYPE,2,QuadraticMeshComponentNumber,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(DependentField,CMFE_FIELD_U_VARIABLE_TYPE,3,QuadraticMeshComponentNumber,Err)
  !CALL cmfe_Field_ComponentMeshComponentSet(DependentField,CMFE_FIELD_U_VARIABLE_TYPE,4,QuadraticMeshComponentNumber,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(DependentField,CMFE_FIELD_U_VARIABLE_TYPE,4,LinearMeshComponentNumber,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(DependentField,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,1,QuadraticMeshComponentNumber,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(DependentField,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,2,QuadraticMeshComponentNumber,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(DependentField,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,3,QuadraticMeshComponentNumber,Err)
  !CALL cmfe_Field_ComponentMeshComponentSet(DependentField,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,4,QuadraticMeshComponentNumber,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(DependentField,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,4,LinearMeshComponentNumber,Err)  
  !                                                 >>>> END 
  CALL cmfe_Field_VariableLabelSet(DependentField,CMFE_FIELD_U_VARIABLE_TYPE,"Dependent",Err)

  !CALL cmfe_Field_NumberOfComponentsSet(DependentField,CMFE_FIELD_V_VARIABLE_TYPE,3,Err)
  !IF(UsePressureBasis) THEN
    ! Set the pressure to be nodally based and use the second mesh component if required
    !CALL cmfe_Field_ComponentInterpolationSet(DependentField,CMFE_FIELD_U_VARIABLE_TYPE,4,CMFE_FIELD_NODE_BASED_INTERPOLATION,Err)
    !CALL cmfe_Field_ComponentInterpolationSet(DependentField,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,4, &
      !& CMFE_FIELD_NODE_BASED_INTERPOLATION,Err)
    !CALL cmfe_Field_ComponentMeshComponentSet(DependentField,CMFE_FIELD_U_VARIABLE_TYPE,4,2,Err)
    !CALL cmfe_Field_ComponentMeshComponentSet(DependentField,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,4,2,Err)
  !END IF
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
  CALL cmfe_Field_ComponentValuesInitialise(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,MPc*C(1),Err)
  CALL cmfe_Field_ComponentValuesInitialise(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,2,MPc*C(2),Err)
  CALL cmfe_Field_ComponentValuesInitialise(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,3,MPc*C(3),Err)
  CALL cmfe_Field_ComponentValuesInitialise(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,4,MPc*C(4),Err)
  CALL cmfe_Field_ComponentValuesInitialise(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,5,MPc*C(5),Err)
  CALL cmfe_Field_ComponentValuesInitialise(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,6,MPc*C(6),Err)
  CALL cmfe_Field_ComponentValuesInitialise(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,7,MPc*C(7),Err)
  CALL cmfe_Field_ComponentValuesInitialise(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,8,MPc*C(8),Err)
  CALL cmfe_Field_ComponentValuesInitialise(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,9,C(9),Err)
  CALL cmfe_Field_ComponentValuesInitialise(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,10,C(10),Err)
  CALL cmfe_Field_ComponentValuesInitialise(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,11,C(11),Err)

  !change the material in the fat/skin to passive
  do elem_idx=13,39
    CALL cmfe_Decomposition_ElementDomainGet(Decomposition,elem_idx,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN
      CALL cmfe_Field_ParameterSetUpdateElement(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
        & elem_idx,10,0.0_CMISSRP,Err)
    ENDIF
  enddo

  ! Create the equations set equations
  CALL cmfe_Equations_Initialise(Equations,Err)
  CALL cmfe_EquationsSet_EquationsCreateStart(EquationsSet,Equations,Err)
  CALL cmfe_Equations_SparsityTypeSet(Equations,CMFE_EQUATIONS_SPARSE_MATRICES,Err)
  CALL cmfe_Equations_OutputTypeSet(Equations,CMFE_EQUATIONS_NO_OUTPUT,Err)
  CALL cmfe_EquationsSet_EquationsCreateFinish(EquationsSet,Err)

  ! Initialise dependent field from undeformed geometry and displacement bcs and set hydrostatic pressure
  CALL cmfe_Field_ParametersToFieldParametersComponentCopy(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 1,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,Err)
  CALL cmfe_Field_ParametersToFieldParametersComponentCopy(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 2,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,2,Err)
  CALL cmfe_Field_ParametersToFieldParametersComponentCopy(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 3,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,3,Err)
  CALL cmfe_Field_ComponentValuesInitialise(DependentField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,4,0.0_CMISSRP,&
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
  !
  ! From TA Example                                >>>> START 
  ! set the main control loop (time loop type)
  !CALL cmfe_ControlLoop_Initialise(ControlLoopMain,Err)
  !CALL cmfe_Problem_ControlLoopGet(Problem,CMFE_CONTROL_LOOP_NODE,ControlLoopMain,Err)
  !CALL cmfe_ControlLoop_LabelSet(ControlLoopMain,'MAIN_TIME_LOOP',Err)
  !CALL cmfe_ControlLoop_TypeSet(ControlLoopMain,CMFE_PROBLEM_CONTROL_TIME_LOOP_TYPE,Err)
  !Loop in time for STIM_STOP with the Stimulus applied.
  !CALL cmfe_ControlLoop_TimesSet(ControlLoopMain,0.0_CMISSRP,ELASTICITY_TIME_STEP,ELASTICITY_TIME_STEP,Err)
  !CALL cmfe_ControlLoop_TimeOutputSet(ControlLoopMain,2,Err)
  !CALL cmfe_ControlLoop_OutputTypeSet(ControlLoopMain,CMFE_CONTROL_LOOP_TIMING_OUTPUT,Err)  
  !                                                >>>> END
  ! control loop simple type <
  !CALL cmfe_ControlLoop_Initialise(ControlLoop,Err)
  !CALL cmfe_Problem_ControlLoopGet(Problem,CMFE_CONTROL_LOOP_NODE,ControlLoop,Err)
  !CALL cmfe_ControlLoop_TypeSet(ControlLoop,CMFE_PROBLEM_CONTROL_LOAD_INCREMENT_LOOP_TYPE,Err)
  !CALL cmfe_ControlLoop_MaximumIterationsSet(ControlLoop,10,Err)
  !CALL cmfe_ControlLoop_LoadOutputSet(ControlLoop,1,Err)
  ! >
  CALL cmfe_Problem_ControlLoopCreateFinish(Problem,Err)

  !Create the problem solvers
  CALL cmfe_Solver_Initialise(Solver,Err)
  CALL cmfe_Solver_Initialise(LinearSolver,Err)
  CALL cmfe_Problem_SolversCreateStart(Problem,Err)
  CALL cmfe_Problem_SolverGet(Problem,CMFE_CONTROL_LOOP_NODE,1,Solver,Err)
  CALL cmfe_Solver_OutputTypeSet(Solver,CMFE_SOLVER_PROGRESS_OUTPUT,Err)
  CALL cmfe_Solver_NewtonJacobianCalculationTypeSet(Solver,CMFE_SOLVER_NEWTON_JACOBIAN_EQUATIONS_CALCULATED,Err)
  CALL cmfe_Solver_NewtonSolutionToleranceSet(Solver,1.E-6_CMISSRP,Err)
  CALL cmfe_Solver_NewtonAbsoluteToleranceSet(Solver,1.E-6_CMISSRP,Err)
  CALL cmfe_Solver_NewtonMaximumIterationsSet(Solver,100,Err)
  CALL cmfe_Solver_NewtonLinearSolverGet(Solver,LinearSolver,Err)
  CALL cmfe_Solver_LinearTypeSet(LinearSolver,CMFE_SOLVER_LINEAR_DIRECT_SOLVE_TYPE,Err)
  !CALL cmfe_Solver_LinearTypeSet(LinearSolver,CMFE_SOLVER_LINEAR_ITERATIVE_SOLVE_TYPE,Err)  
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

!Fix TOP_NODES in all directions
VALUE = 0.0_CMISSRP
  do i=1,21
    NodeNumber=TOP_NODES(i)
    CALL cmfe_Decomposition_NodeDomainGet(Decomposition,NodeNumber,1,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN
      CALL cmfe_Field_ParameterSetGetNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,NodeNumber, &
       & 1,VALUE,Err)
      CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,NodeNumber,1, &
       & CMFE_BOUNDARY_CONDITION_FIXED,VALUE,Err)
      CALL cmfe_Field_ParameterSetGetNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,NodeNumber, & 
       & 2,VALUE,Err)
      CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,NodeNumber,2, &
       & CMFE_BOUNDARY_CONDITION_FIXED,VALUE,Err)
      CALL cmfe_Field_ParameterSetGetNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,NodeNumber, &
       & 3,VALUE,Err)
      CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,NodeNumber,3, &
       & CMFE_BOUNDARY_CONDITION_FIXED,VALUE,Err)
    ENDIF
  enddo

  !Fix BOTTOM_NODES in all directions
  do i=1,21
    NodeNumber=BOTTOM_NODES(i)
    CALL cmfe_Decomposition_NodeDomainGet(Decomposition,NodeNumber,1,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN
      CALL cmfe_Field_ParameterSetGetNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,NodeNumber, &
       & 1,VALUE,Err)
      CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,NodeNumber,1, &
       & CMFE_BOUNDARY_CONDITION_FIXED,VALUE,Err)
      CALL cmfe_Field_ParameterSetGetNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,NodeNumber, & 
       & 2,VALUE,Err)
      CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,NodeNumber,2, &
       & CMFE_BOUNDARY_CONDITION_FIXED,VALUE,Err)
      CALL cmfe_Field_ParameterSetGetNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,NodeNumber, &
       & 3,VALUE,Err)
      CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,NodeNumber,3, &
       & CMFE_BOUNDARY_CONDITION_FIXED,VALUE,Err)
    ENDIF
  enddo

  !Fix SIDE_NODES in all directions
  do i=1,52
    NodeNumber=SIDE_NODES(i)
    CALL cmfe_Decomposition_NodeDomainGet(Decomposition,NodeNumber,1,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN
      CALL cmfe_Field_ParameterSetGetNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,NodeNumber, &
       & 1,VALUE,Err)
      CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,NodeNumber,1, &
       & CMFE_BOUNDARY_CONDITION_FIXED,VALUE,Err)
      CALL cmfe_Field_ParameterSetGetNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,NodeNumber, & 
       & 2,VALUE,Err)
      CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,NodeNumber,2, &
       & CMFE_BOUNDARY_CONDITION_FIXED,VALUE,Err)
      CALL cmfe_Field_ParameterSetGetNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,NodeNumber, &
       & 3,VALUE,Err)
      CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,NodeNumber,3, &
       & CMFE_BOUNDARY_CONDITION_FIXED,VALUE,Err)
    ENDIF
  enddo
write(*,*) "applied BCs"
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

  CALL cmfe_SolverEquations_BoundaryConditionsCreateFinish(SolverEquations,Err)

!
! Solving the problem---------------------------------------------------------------------------------------------------------START
!
  !Output solution
  CALL cmfe_Fields_Initialise(Fields,Err)
  CALL cmfe_Fields_Create(Region,Fields,Err)  
  CALL cmfe_Fields_NodesExport(Fields,"ActiveStrain_TransIso","FORTRAN",Err)
  CALL cmfe_Fields_ElementsExport(Fields,"ActiveStrain_TransIso","FORTRAN",Err)
  CALL cmfe_Fields_Finalise(Fields,Err)

  ! Read in activation alpha at time t; no check currently made that TIMESTEP = number of rows in text file
  !OPEN (UNIT=3, FILE='activation_alpha.txt', STATUS='OLD', ACTION='read')
  !DO i=1,TIMESTEPS
    !READ(3,*) alpha(i)
  !ENDDO
  !CLOSE(3)
  
  ! Initialise alpha(i=t) 
  !CALL cmfe_Field_ComponentValuesInitialise(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,9, &
  !  & 0.0_CMISSRP,Err)
 
  ! Output solution
  CALL cmfe_Fields_Initialise(Fields,Err)
  CALL cmfe_Fields_Create(Region,Fields,Err)  

  ! Solve the mechanical problem fpr this time step - DUMMY SOLVE - required to initialise all variables 
  ! (ideally should not be required)
  !CALL cmfe_Problem_Solve(Problem,Err)

  ! Mechanical BC increment (if concentric contraction is used)
  BCLOAD = 0.0_CMISSRP
  time = 0.0_CMISSRP
  alpha_t = 0.0_CMISSRP
  i = 0
  ! Loop over Time -----------------------------------------------------------------------------------------------------------START
  DO WHILE(time <= TIME_STOP)
    WRITE(*,*) "-------------------------------------------------------"
    WRITE(*,*) "Percentage completed: ", 100.0_CMISSRP*time/TIME_STOP
    
    ! Solve the mechanical problem for this time step - DUMMY SOLVE - required to initialise all variables 
    ! (ideally should not be required)
    !CALL cmfe_Problem_Solve(Problem,Err)

    ! Below are different options to update alpha at each time step ----------------------------------------------------------START
      ! (1) Entire muscle (working)
      ! (2) Element based (working)
      ! (3) Node based
      ! (4) Gauss-point based 

      ! Set spatial update case - WARNING: cmfe_Field_ComponentInterpolationSet for C(9) needs to be correspondingly set!  
      alpha_SU = 1

      SELECT CASE (alpha_SU)
        CASE(1) ! Entire muscle (working)
          CALL cmfe_Field_ComponentValuesInitialise(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,9, &
            & alpha_t,Err)
          write(*,*) "Activation, MaxStress, Activation*MaxStress: ", alpha_t, C(11), alpha_t*C(11)
        !CASE(2) ! Element based (working)
          ! loop over all elements
          !DO elem_idx=1,NumberOfElementsFE
            !CALL cmfe_Field_ParameterSetUpdateElement(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,&
              !& elem_idx,9,alpha_t,Err)
          !END DO        
        
        !CASE(3) ! Node based
          ! loop over all nodes
          !DO node_idx=1,TotalNumberOfNodes 
            !CALL cmfe_Field_ParameterSetUpdateNode(MaterialField,& 
              !& CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1, 1, node_idx,9,alpha_t,Err)      
          !END DO       

        !CASE(4) !Gauss-point based
          ! loop over all elements
          !DO elem_idx=1,NumberOfElementsFE
            ! loop over all gauss points
            !DO gauss_idx = 1, NumberOfGaussXi
              !CALL cmfe_Field_ParameterSetUpdateGaussPoint(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,&
                !& gauss_idx, elem_idx, 9, alpha_t, Err)
            !END DO 
          !END DO       

        CASE DEFAULT ! Entire muscle (working)
           CALL cmfe_Field_ComponentValuesInitialise(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,9, &
            & alpha_t,Err)     
      END SELECT! Spatial alpha update------------------------------------------------------------------------------------------END
    
    ! time loop
    !CALL cmfe_ControlLoop_TimesSet(ControlLoopMain,time,time+PERIOD,ELASTICITY_TIME_STEP,Err)
    ! Solve the mechanical problem fpr this time step
    CALL cmfe_Problem_Solve(Problem,Err)
    time = time + PERIOD
  
    ! Set how alpha evolves over time (may choose to keep it constant, scale it, etc.)
    !alpha_t = exp(alpha(i))-1.0_CMISSRP
    i = i+1
    alpha_t = alpha_t + alpha_inc
    ! update the mechanical boundary condition -------------------------------------------------------------------------------START
      !DO node_idx=1,SIZE(BOTTOM_NODES,1)
        !NodeNumber=BOTTOM_NODES(node_idx)
        !CALL cmfe_Decomposition_NodeDomainGet(Decomposition,NodeNumber,1,NodeDomain,Err)
        !IF(NodeDomain==ComputationalNodeNumber) THEN
          !CALL cmfe_Field_ParameterSetAddNode(DependentField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, & 
            !& 1,1,NodeNumber,3, BCLOAD,Err)
        !ENDIF
      !ENDDO ! mechanical BC update ---------------------------------------------------------------------------------------------END
 
    ! Output the exnode file at current time step 
    ! NOTE: if changing filename remember to change AXY!
    IF (i.LE.9) THEN
      WRITE(filename, "(A22,I1)") "ActiveStrain_TransIso_", i
      filename=trim(filename)
    ELSEIF (i.LE.99) THEN
      WRITE(filename, "(A22,I2)") "ActiveStrain_TransIso_", i
      filename=trim(filename)
    ELSE  
      WRITE(filename, "(A22,I3)") "ActiveStrain_TransIso_", i
      filename=trim(filename)
    ENDIF
    CALL cmfe_Fields_NodesExport(Fields,filename,"FORTRAN",Err)
  
  END DO
  ! Loop over Time -------------------------------------------------------------------------------------------------------------END  
  
  !Output solution
  CALL cmfe_Fields_NodesExport(Fields,"ActiveStrain_TransIso","FORTRAN",Err)
  CALL cmfe_Fields_ElementsExport(Fields,"ActiveStrain_TransIso","FORTRAN",Err)
  CALL cmfe_Fields_Finalise(Fields,Err)

!
! ---------------------------------------------------------------------------------------------------------------------------------
!

  CALL cmfe_Finalise(Err)

  WRITE(*,'(A)') "Program successfully completed."

  STOP

END PROGRAM LARGEUNIAXIALEXTENSIONEXAMPLE



