# Example application target
file(GLOB exe_sources *.*90)
add_executable(
  "${PROJECT_NAME}-exe"
  ${exe_sources}
)
# executable name should be project name
set_target_properties(
  "${PROJECT_NAME}-exe"
  PROPERTIES
  OUTPUT_NAME "${PROJECT_NAME}"
)
# link against the library in src/
target_link_libraries(
  "${PROJECT_NAME}-exe"
  PRIVATE
  "${PROJECT_NAME}-lib"
)
# make the app instalable
install(
  TARGETS
  "${PROJECT_NAME}-exe"
  RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}"
)
