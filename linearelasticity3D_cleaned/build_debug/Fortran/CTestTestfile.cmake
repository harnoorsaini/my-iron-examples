# CMake generated Testfile for 
# Source directory: /usr/local/opencmiss/examples/linearelasticity3D_cleaned/Fortran
# Build directory: /usr/local/opencmiss/examples/linearelasticity3D_cleaned/build_debug/Fortran
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(Uniaxial_3D "/usr/local/opencmiss/examples/linearelasticity3D_cleaned/build_debug/Fortran/Uniaxial_3D")
set_tests_properties(Uniaxial_3D PROPERTIES  ENVIRONMENT "LD_LIBRARY_PATH=/usr/local/opencmiss/iron/install/x86_64_linux/gnu-5.4-F5.4/openmpi_release/release/bin:")
