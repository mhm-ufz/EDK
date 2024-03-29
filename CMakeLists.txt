cmake_minimum_required(VERSION 3.14)
# get version and date from files (version.txt and version_date.txt)
include(cmake/version.cmake)
get_version(EDK_VER EDK_VER_DEV EDK_DATE)

# define the project
project(
  edk
  VERSION ${EDK_VER}
  LANGUAGES Fortran C
)

# add the library
add_subdirectory(src)
set(LIB_NAME "${PROJECT_NAME}-lib")

# add the executable
add_subdirectory(app)

# add full version and date to pre-processor flags (qoutes need in before hand)
target_compile_definitions(${LIB_NAME} PRIVATE
  EDKVERSION='${EDK_VER_DEV}'
  EDKDATE='${EDK_DATE}'
)

# prepare coverage
if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU" AND CMAKE_WITH_COVERAGE)
  message(STATUS "${PROJECT_NAME}: coverage enabled")
  include(cmake/CodeCoverage.cmake)
  append_coverage_compiler_flags_to_target(${LIB_NAME})
  append_coverage_compiler_flags_to_target(${PROJECT_NAME}-exe)
  SETUP_TARGET_FOR_COVERAGE_LCOV(
    NAME coverage
    EXECUTABLE ../scripts/run_cmake_coverage.sh
    DEPENDENCIES ${LIB_NAME} ${PROJECT_NAME}-exe
    EXCLUDE "src/tests/*" "${CMAKE_CURRENT_BINARY_DIR}/*"
    LCOV_ARGS --no-external
    GENHTML_ARGS -t "${PROJECT_NAME} coverage" --html-prolog ../doc/html_files/cov_header.prolog
  )
endif()
