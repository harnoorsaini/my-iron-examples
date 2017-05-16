# CMake generated Testfile for 
# Source directory: /usr/local/opencmiss/examples/input_file_example/Fortran
# Build directory: /usr/local/opencmiss/examples/input_file_example/build_debug/Fortran
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(Iron_GE "/usr/local/opencmiss/examples/input_file_example/build_debug/Fortran/iron_ge")
set_tests_properties(Iron_GE PROPERTIES  ENVIRONMENT "LD_LIBRARY_PATH=/usr/local/opencmiss/iron/install/x86_64_linux/gnu-5.4-F5.4/openmpi_release/release/bin:")
