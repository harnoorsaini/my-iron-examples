# CMake generated Testfile for 
# Source directory: /usr/local/opencmiss/examples/user_mesh/Fortran
# Build directory: /usr/local/opencmiss/examples/user_mesh/build_debug/Fortran
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(MoreComplexMeshExample "/usr/local/opencmiss/examples/user_mesh/build_debug/Fortran/MoreComplexMeshExample")
set_tests_properties(MoreComplexMeshExample PROPERTIES  ENVIRONMENT "LD_LIBRARY_PATH=/usr/local/opencmiss/iron/install/x86_64_linux/gnu-5.4-F5.4/openmpi_release/release/bin:")
