# CMake generated Testfile for 
# Source directory: /usr/local/opencmiss/examples/simple_muscle_harry/Fortran
# Build directory: /usr/local/opencmiss/examples/simple_muscle_harry/build_debug/Fortran
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(Simple_ActiveStrainHarry "/usr/local/opencmiss/examples/simple_muscle_harry/build_debug/Fortran/Simple_ActiveStrainHarry")
set_tests_properties(Simple_ActiveStrainHarry PROPERTIES  ENVIRONMENT "LD_LIBRARY_PATH=/usr/local/opencmiss/iron/install/x86_64_linux/gnu-5.4-F5.4/openmpi_release/release/bin:")
