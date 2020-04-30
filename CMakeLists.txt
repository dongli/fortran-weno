cmake_minimum_required(VERSION 3.0)

project(fortran-weno LANGUAGES Fortran)

set(source_files
  src/weno_interp_mod.F90
  src/weno.F90
)

add_library(weno ${source_files})
