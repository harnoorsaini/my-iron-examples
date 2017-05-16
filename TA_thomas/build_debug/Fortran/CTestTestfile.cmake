# CMake generated Testfile for 
# Source directory: /usr/local/opencmiss/examples/TA_thomas/Fortran
# Build directory: /usr/local/opencmiss/examples/TA_thomas/build_debug/Fortran
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(TAexample "/usr/local/opencmiss/examples/TA_thomas/build_debug/Fortran/TAexample")
set_tests_properties(TAexample PROPERTIES  ENVIRONMENT "LD_LIBRARY_PATH=/usr/local/opencmiss/iron/install/x86_64_linux/gnu-5.4-F5.4/openmpi_release/release/bin:")
