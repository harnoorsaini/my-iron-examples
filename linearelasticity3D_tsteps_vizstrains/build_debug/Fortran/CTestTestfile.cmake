# CMake generated Testfile for 
# Source directory: /usr/local/opencmiss/examples/linearelasticity3D_tsteps_vizstrains/Fortran
# Build directory: /usr/local/opencmiss/examples/linearelasticity3D_tsteps_vizstrains/build_debug/Fortran
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(linelastic2Dtstep_strain "/usr/local/opencmiss/examples/linearelasticity3D_tsteps_vizstrains/build_debug/Fortran/linelastic2Dtstep_strain")
set_tests_properties(linelastic2Dtstep_strain PROPERTIES  ENVIRONMENT "LD_LIBRARY_PATH=/usr/local/opencmiss/iron/install/x86_64_linux/gnu-5.4-F5.4/openmpi_release/release/bin:")
