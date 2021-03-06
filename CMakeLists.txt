cmake_minimum_required (VERSION 2.8)
project (YIP)
enable_language (Fortran)

# Default to Release
if (NOT CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build, options are: None Debug Release." FORCE)
endif (NOT CMAKE_BUILD_TYPE)

# Require MPI
find_package(MPI REQUIRED)
set( CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})

# Setup Compile Flags
set (CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} -i4 -r8")
set (CMAKE_Fortran_FLAGS_DEBUG   "${CMAKE_Fortran_FLAGS_DEBUG} -g -traceback -check pointer  -check bounds -i4 -r8")

# Build Library
add_library (fft STATIC "fft.f")
add_library (kdtree STATIC "kdtree.f90")

# Build Executable
file(GLOB MOD Modules/*.F)
file(GLOB SRC Subroutines/*.F)
set_source_files_properties(${MOD} ${SRC} COMPILE_FLAGS "-assume byterecl")

add_executable (lesmpi.a ${SRC} ${MOD})
target_link_libraries (lesmpi.a fft kdtree)

install ( TARGETS lesmpi.a RUNTIME DESTINATION "bin")
