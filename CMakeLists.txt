cmake_minimum_required(VERSION 3.14)
# get version and date from files (version.txt and version_date.txt)
include(cmake/cmake-modules/version.cmake)
get_version(EDK_VER EDK_VER_DEV EDK_DATE)

# define the project ###########################################################################
project(                                                                                       #
  edk                                                                                          #
  VERSION ${EDK_VER}                                                                           #
  LANGUAGES Fortran                                                                            #
)                                                                                              #
################################################################################################

# add the library
add_subdirectory(./src)
set(LIB_NAME "${PROJECT_NAME}-lib")

# (optional) add the executable
option(BUILD_${PROJECT_NAME}_APP "Build project with app." ON)
if(EXISTS ./app AND BUILD_${PROJECT_NAME}_APP)
  message(STATUS "${PROJECT_NAME}: found app")
  add_subdirectory(./app)
endif()

# add full version and date to pre-processor flags (qoutes need in before hand)
target_compile_definitions(${LIB_NAME} PRIVATE
  EDKVERSION='${EDK_VER_DEV}'
  EDKDATE='${EDK_DATE}'
)

# automatically enables testing (OFF by default)
option(BUILD_TESTING "Build with (pfUnit) tests." ${CMAKE_WITH_COVERAGE})
include(CTest)

# add test directories
if(BUILD_TESTING)
  message(STATUS "${PROJECT_NAME}: build tests")
  add_subdirectory(./src/tests)
endif()

# prepare coverage
if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU" AND CMAKE_WITH_COVERAGE AND BUILD_TESTING)
  message(STATUS "${PROJECT_NAME}: coverage enabled")
  file(GLOB testfiles ./src/tests/*.pf)
  if(NOT testfiles)
    message(STATUS "${PROJECT_NAME}: no tests found for coverage")
    set(_add_lcov_args -i)
  endif()
  include(cmake/cmake-modules/CodeCoverage.cmake)
  target_compile_options(${LIB_NAME} PRIVATE -g --coverage)
  target_link_libraries(${LIB_NAME} PRIVATE gcov)
  SETUP_TARGET_FOR_COVERAGE_LCOV(
    NAME coverage
    EXECUTABLE ctest
    DEPENDENCIES ${LIB_NAME}
    EXCLUDE "app/*" "src/tests/*" "${CMAKE_CURRENT_BINARY_DIR}/*"
    LCOV_ARGS --no-external ${_add_lcov_args}
    GENHTML_ARGS -t "${PROJECT_NAME} coverage" --html-prolog ../doc/html_files/cov_header.prolog
  )
endif()
