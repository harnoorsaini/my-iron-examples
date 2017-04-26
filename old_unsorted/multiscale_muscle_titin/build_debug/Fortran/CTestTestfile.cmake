# CMake generated Testfile for 
# Source directory: /usr/local/opencmiss/examples/multiscale_muscle_titin/Fortran
# Build directory: /usr/local/opencmiss/examples/multiscale_muscle_titin/build_debug/Fortran
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(Laplace_Fortran "/usr/local/opencmiss/examples/multiscale_muscle_titin/build_debug/Fortran/laplace_fortran")
set_tests_properties(Laplace_Fortran PROPERTIES  ENVIRONMENT "LD_LIBRARY_PATH=/usr/local/opencmiss/iron/install/x86_64_linux/gnu-5.4-F5.4/openmpi_release/debug/bin:")
