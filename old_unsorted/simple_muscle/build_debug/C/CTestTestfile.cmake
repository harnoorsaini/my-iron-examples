# CMake generated Testfile for 
# Source directory: /usr/local/opencmiss/examples/simple_muscle/C
# Build directory: /usr/local/opencmiss/examples/simple_muscle/build_debug/C
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(Laplace_C "laplace_c")
set_tests_properties(Laplace_C PROPERTIES  ENVIRONMENT "LD_LIBRARY_PATH=/usr/local/opencmiss/iron/install/x86_64_linux/gnu-5.4-F5.4/openmpi_release/debug/bin:")
