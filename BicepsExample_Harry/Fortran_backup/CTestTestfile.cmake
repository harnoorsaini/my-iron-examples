# CMake generated Testfile for 
# Source directory: /usr/local/opencmiss/examples/BicepsExample_Harry/Fortran
# Build directory: /usr/local/opencmiss/examples/BicepsExample_Harry/build_debug/Fortran
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(BicepsExample_Harry "/usr/local/opencmiss/examples/BicepsExample_Harry/build_debug/Fortran/BicepsExample_Harry")
set_tests_properties(BicepsExample_Harry PROPERTIES  ENVIRONMENT "LD_LIBRARY_PATH=/usr/local/opencmiss/iron/install/x86_64_linux/gnu-5.4-F5.4/openmpi_release/debug/bin:")
