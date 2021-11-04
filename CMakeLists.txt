cmake_minimum_required(VERSION 3.0)

project(fortran-weno LANGUAGES Fortran)

if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffree-line-length-none -ffpe-trap=invalid,zero,overflow")
  if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O0 -fcheck=all -g -fbacktrace")
  else ()
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Ofast -lgomp -fopenmp")
  endif ()
elseif (CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
  if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -no-wrap-margin -O0 -check bounds -g -traceback")
  else ()
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -no-wrap-margin -O3 -qopenmp")
  endif ()
endif ()
set(CMAKE_Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}")

set(source_files
  src/math_mod.F90
  src/poly_utils_mod.F90
  src/smooth_indicators_mod.F90
  src/weno_tensor_product_mod.F90
  src/weno_rhombus_mod.F90
  src/poly_tensor_product_mod.F90
  src/weno.F90
)

add_library(weno ${source_files})

get_directory_property(parent_dir PARENT_DIRECTORY)
if (NOT parent_dir)
  add_executable(test.exe src/test.f90)
  target_link_libraries(test.exe weno)
endif ()
