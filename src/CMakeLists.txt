cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

project(FFTsurf)

enable_language(Fortran)


# this is where we will place the Fortran module files
set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/modules)


add_subdirectory(intersection) # here we use intersection/CMakeLists.txt


add_executable(
  fftsurf.out
  fftsurf.f90
  modules/mod_obb.f90
  modules/mod_chebyshev.f90
  modules/mod_separation.f90
  modules/mod_diffgeom.f90
  modules/mod_tolerances.f90
  modules/mod_errors_intersection.f90
  modules/mod_types_dcel.f90
  modules/mod_linearprogramming.f90
  modules/mod_types_intersection.f90
  modules/mod_math.f90
  modules/mod_util.f90
)

target_link_libraries(
  fftsurf.out ${PROJECT_SOURCE_DIR}/dfftpack/libdfftpack.a
  )
