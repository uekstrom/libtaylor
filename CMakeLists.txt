cmake_minimum_required(VERSION 3.5 FATAL_ERROR)
project(taylor
  LANGUAGES CXX
  VERSION 1.0.0
  )

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

include(GNUInstallDirs)
set(INSTALL_INCLUDEDIR ${CMAKE_INSTALL_INCLUDEDIR} CACHE PATH "Installation directory for header files")
if(WIN32 AND NOT CYGWIN)
  set(DEF_INSTALL_CMAKEDIR CMake)
else()
  set(DEF_INSTALL_CMAKEDIR share/cmake/${PROJECT_NAME})
endif()
set(INSTALL_CMAKEDIR ${DEF_INSTALL_CMAKEDIR} CACHE PATH "Installation directory for CMake files")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_BINDIR})

enable_testing()

add_subdirectory(tests)

# <<< Create a target >>> Install and export targets >>>

add_library(taylor INTERFACE)

target_compile_definitions(taylor
  INTERFACE
    $<$<BOOL:${WIN32}>:_USE_MATH_DEFINES>
  )

target_include_directories(taylor
  INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:${INSTALL_INCLUDEDIR}>
  )

# <<< Install header files >>>

install(
  FILES
    ctaylor.hpp
    ctaylor_math.hpp
    polymul.hpp
    taylor.hpp
    taylor_math.hpp
    tmath.hpp
  DESTINATION
    ${INSTALL_INCLUDEDIR}/taylor
  )

# <<< Export and install target >>>

install(
  TARGETS
    taylor
  EXPORT
    taylorTargets
  )

install(
  EXPORT
    taylorTargets
  NAMESPACE
    "taylor::"
  DESTINATION
    ${INSTALL_CMAKEDIR}
  COMPONENT
    dev
  )

include(CMakePackageConfigHelpers)
write_basic_package_version_file(
  ${CMAKE_CURRENT_BINARY_DIR}/taylorConfigVersion.cmake
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY SameMajorVersion
  )

configure_package_config_file(
  ${PROJECT_SOURCE_DIR}/cmake/taylorConfig.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/taylorConfig.cmake
  INSTALL_DESTINATION ${INSTALL_CMAKEDIR}
  )

install(
  FILES
    ${CMAKE_CURRENT_BINARY_DIR}/taylorConfig.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/taylorConfigVersion.cmake
  DESTINATION
    ${INSTALL_CMAKEDIR}
  )

set_property(DIRECTORY PROPERTY EP_BASE ${CMAKE_BINARY_DIR}/subprojects)

# we need ExternalProject functionality
include(ExternalProject)

# fetch pybind11 sources
ExternalProject_Add(pybind11_external
  GIT_REPOSITORY
    https://github.com/pybind/pybind11.git
  GIT_TAG
    v2.2.3
  CMAKE_ARGS
    -DCMAKE_INSTALL_PREFIX=${CMAKE_CURRENT_BINARY_DIR}/stage
    -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
    -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
    -DCMAKE_CXX_STANDARD=${CMAKE_CXX_STANDARD}
    -DPYBIND11_TEST=OFF
  LOG_DOWNLOAD 1
  LOG_UPDATE 1
  LOG_CONFIGURE 1
  LOG_BUILD 1
  LOG_INSTALL 1
  )

# process our source code as if it was an external project
ExternalProject_Add(pytaylor
  DEPENDS pybind11_external
  SOURCE_DIR
    ${CMAKE_CURRENT_SOURCE_DIR}/pytaylor
  BINARY_DIR
    ${CMAKE_CURRENT_BINARY_DIR}/pytaylor
  CMAKE_ARGS
    -Dpybind11_DIR=${CMAKE_CURRENT_BINARY_DIR}/stage/share/cmake/pybind11
    -DCMAKE_INSTALL_PREFIX=${CMAKE_CURRENT_BINARY_DIR}/stage
    -DLIBTAYLOR_INCLUDE_DIR=${CMAKE_CURRENT_SOURCE_DIR}
    -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
    -DCMAKE_CXX_STANDARD=${CMAKE_CXX_STANDARD}
  )

# define test
add_test(
  NAME
    python_test
  COMMAND
    ${CMAKE_COMMAND} -E env PYTAYLOR_MODULE_PATH=${CMAKE_CURRENT_BINARY_DIR}
    ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/pytaylor/test_taylor_import.py
    #${PYTHON_EXECUTABLE} -m pytest -rws -v --capture=sys ${CMAKE_CURRENT_SOURCE_DIR}/pytaylor/test.py
  WORKING_DIRECTORY
    ${CMAKE_CURRENT_BINARY_DIR}
  )
