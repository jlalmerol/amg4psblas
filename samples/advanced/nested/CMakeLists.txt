cmake_minimum_required(VERSION 3.15)
project(AMGNestedExamples Fortran)

set(AMG4PSBLAS_INSTALL_DIR "" CACHE PATH "Path to AMG installation")
set(PSBLAS_INSTALL_DIR "" CACHE PATH "Path to PSBLAS installation")

if(NOT AMG4PSBLAS_INSTALL_DIR)
  message(FATAL_ERROR "AMG4PSBLAS_INSTALL_DIR must be set.")
endif()

if(NOT PSBLAS_INSTALL_DIR)
  message(FATAL_ERROR "PSBLAS_INSTALL_DIR must be set.")
endif()

find_package(psblas REQUIRED PATHS ${PSBLAS_INSTALL_DIR})
find_package(amg4psblas REQUIRED PATHS ${AMG4PSBLAS_INSTALL_DIR})
find_package(MPI REQUIRED COMPONENTS Fortran)

include_directories(
  "."
  "${AMG4PSBLAS_INSTALL_DIR}/include"
  "${AMG4PSBLAS_INSTALL_DIR}/modules"
  "${PSBLAS_INSTALL_DIR}/include"
  "${PSBLAS_INSTALL_DIR}/modules"
)

set(CMAKE_Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/modules")
set(EXEDIR "${CMAKE_CURRENT_BINARY_DIR}/runs")
file(MAKE_DIRECTORY "${EXEDIR}")

set(AMG_LIBS amg4psblas::amgprec)
set(PSBLAS_LIBS psblas::util psblas::linsolve psblas::prec psblas::base)

add_library(amg_d_nested_block_prec STATIC amg_d_nested_block_prec_mod.F90)
target_link_libraries(amg_d_nested_block_prec PUBLIC
  ${AMG_LIBS}
  ${PSBLAS_LIBS}
  MPI::MPI_Fortran
  stdc++
)

function(create_nested_executable target source)
  add_executable(${target} ${source})
  target_link_libraries(${target} PUBLIC amg_d_nested_block_prec)
  set_target_properties(${target} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${EXEDIR})
endfunction()

create_nested_executable(amg_d_nest_stokes_file_test amg_d_nest_stokes_file_test.F90)
create_nested_executable(amg_d_nest_optimal_control_test amg_d_nest_optimal_control_test.F90)

install(TARGETS
  amg_d_nest_stokes_file_test
  amg_d_nest_optimal_control_test
  DESTINATION bin
)
