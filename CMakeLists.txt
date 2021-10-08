cmake_minimum_required(VERSION 3.0)

project(fortran-weno LANGUAGES Fortran)

set(source_files
  src/math_mod.F90
  src/poly_utils_mod.F90
  src/weno_types_mod.F90
  src/weno_coef_mod.F90
)

add_library(weno ${source_files})

add_executable(test.exe src/test.f90)
target_link_libraries(test.exe weno)
