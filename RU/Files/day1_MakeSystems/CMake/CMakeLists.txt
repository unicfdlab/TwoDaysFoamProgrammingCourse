#
# Minimum cmake version
#
cmake_minimum_required (VERSION 2.6)

#
# Project name
#
project (test_main)

#
# Project folder
#
include_directories (.)

add_executable (main main.cxx circle.cxx square.cxx fprotos.hxx)

target_link_libraries(main m)

#END_OF_FILE

