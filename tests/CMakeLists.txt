add_executable(unittest_taylor ${CMAKE_CURRENT_SOURCE_DIR}/unittest_taylor.cpp)

target_link_libraries(unittest_taylor
  PRIVATE
    taylor
  )

add_test(
  NAME Taylor
  COMMAND $<TARGET_FILE:unittest_taylor>
  )

# <<< Test CMake module is properly exported

add_test(
  NAME use-taylor_configure
  COMMAND
    ${CMAKE_COMMAND} -H${CMAKE_CURRENT_LIST_DIR}/use_taylor
                     -B${CMAKE_CURRENT_BINARY_DIR}/build_use-taylor
                     -G${CMAKE_GENERATOR}
                     -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
                     -Dtaylor_DIR=${CMAKE_INSTALL_PREFIX}/${INSTALL_CMAKEDIR}
                     -DCMAKE_BUILD_TYPE=$<CONFIGURATION>
  )
set_tests_properties(use-taylor_configure
  PROPERTIES
    FIXTURES_SETUP use-taylor
  )

add_test(
  NAME use-taylor_build
  COMMAND
    ${CMAKE_COMMAND} --build ${CMAKE_CURRENT_BINARY_DIR}/build_use-taylor
                     --config $<CONFIGURATION>
  )
set_tests_properties(use-taylor_build
  PROPERTIES
    FIXTURES_REQUIRED use-taylor
  )

add_test(
  NAME use-taylor_cleanup
  COMMAND
    ${CMAKE_COMMAND} -E remove_directory ${CMAKE_CURRENT_BINARY_DIR}/build_use-taylor
  )
set_tests_properties(use-taylor_cleanup
  PROPERTIES
    FIXTURES_CLEANUP use-taylor
  )

# <<< Test Python bindings >>>
add_test(
  NAME
    pytaylor_import
  COMMAND
   ${CMAKE_COMMAND} -E env PYTAYLOR_MODULE_PATH=${STAGED_INSTALL_PREFIX}/${INSTALL_LIBDIR}
   ${PYTHON_EXECUTABLE} -m pytest -rws -v --capture=sys ${CMAKE_CURRENT_SOURCE_DIR}/test_taylor_import.py
  WORKING_DIRECTORY
    ${CMAKE_CURRENT_BINARY_DIR}
  )
