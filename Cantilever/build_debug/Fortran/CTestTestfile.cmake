# CMake generated Testfile for 
# Source directory: /usr/local/opencmiss/examples/Cantilever/Fortran
# Build directory: /usr/local/opencmiss/examples/Cantilever/build_debug/Fortran
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(Cantilever "/usr/local/opencmiss/examples/Cantilever/build_debug/Fortran/Cantilever")
set_tests_properties(Cantilever PROPERTIES  ENVIRONMENT "LD_LIBRARY_PATH=/usr/local/opencmiss/iron/install/x86_64_linux/gnu-5.4-F5.4/openmpi_release/release/bin:")
