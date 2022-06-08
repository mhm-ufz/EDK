set(LIB_NAME "${PROJECT_NAME}-lib")

# add sources .f90 and .F90 files
file(GLOB sources ./*.*90)
option(BUILD_${PROJECT_NAME}_LIB_SHARED "Build library as shared." OFF)
if(BUILD_${PROJECT_NAME}_LIB_SHARED)
  add_library(${LIB_NAME} SHARED ${sources})
else()
  add_library(${LIB_NAME} STATIC ${sources})
endif()
target_include_directories(${LIB_NAME} PUBLIC ${CMAKE_CURRENT_BINARY_DIR})

# download forces and link it
include(../cmake/cmake-modules/CPM.cmake)
CPMAddPackage("https://git.ufz.de/chs/forces.git@0.3.1")
if(BUILD_${PROJECT_NAME}_LIB_SHARED)
  set_property(TARGET forces PROPERTY POSITION_INDEPENDENT_CODE ON)
endif()
target_link_libraries(${LIB_NAME} PUBLIC forces)

# add all compile options (MPI, OpenMP, Lapack, Coverage)
set(CMAKE_WITH_LAPACK ON)
include(../cmake/cmake-modules/compileoptions.cmake)
if (CMAKE_WITH_OpenMP)
  target_link_libraries(${LIB_NAME} PUBLIC OpenMP::OpenMP_Fortran)
endif()
if (CMAKE_WITH_LAPACK)
  target_link_libraries(${LIB_NAME} PUBLIC ${LAPACK_LIBRARIES})
endif()

# set compiling flags for debug and relese version
if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  target_compile_definitions(${LIB_NAME} PUBLIC "GFORTRAN")
  target_compile_options(${LIB_NAME} PUBLIC
    -ffree-form -ffixed-line-length-132
    $<$<CONFIG:DEBUG>:-Og -g -Wall -Wextra -fimplicit-none -fbacktrace -fcheck=all -ffpe-trap=zero,overflow,underflow -finit-real=snan -pedantic-errors>
    $<$<CONFIG:RELEASE>:-O3>
    $<$<BOOL:${CMAKE_WITH_GPROF}>:-pg>
  )
endif()
if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  target_compile_definitions(${LIB_NAME} PUBLIC "INTEL")
  target_compile_options(${LIB_NAME} PUBLIC
    -nofixed "SHELL:-assume byterecl" "SHELL:-fp-model source" -m64 "SHELL:-assume realloc-lhs"
    $<$<CONFIG:DEBUG>:-g "SHELL:-warn all" "SHELL:-check all" -debug -traceback -fp-stack-check -O0>
    $<$<CONFIG:RELEASE>:-O3 -qoverride-limits>
  )
endif()
if(CMAKE_Fortran_COMPILER_ID MATCHES "NAG")
  target_compile_definitions(${LIB_NAME} PUBLIC "NAG")
  target_compile_options(${LIB_NAME} PUBLIC
    -colour -unsharedf95 -ideclient
    $<$<CONFIG:DEBUG>:-gline -g -nan -O0 -C=all -strict95 -ieee=full>
    $<$<CONFIG:RELEASE>:-O4 -ieee=full>
  )
endif()
