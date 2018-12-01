# define minimum cmake version
cmake_minimum_required(VERSION 3.5 FATAL_ERROR)

# project name and supported language
project(pytaylor
  LANGUAGES CXX
  )

# require C++11
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

#  <<< pybind11 & Python >>>
find_package(pybind11 2.2.3 CONFIG REQUIRED)
message(STATUS "Using pybind11: ${pybind11_INCLUDE_DIR} (version ${pybind11_VERSION} for Py${PYTHON_VERSION_STRING})")
message(STATUS "Using Python ${PYTHON_VERSION_MAJOR}.${PYTHON_VERSION_MINOR}: ${PYTHON_EXECUTABLE}")

#  <<< compiled extension target >>>
pybind11_add_module(pytaylor MODULE THIN_LTO pytaylor.cpp)

if(WIN32 AND (NOT MSVC))
  target_compile_definitions(pytaylor
    PUBLIC
      "_hypot=hypot"
    )
endif()

target_include_directories(pytaylor
  PUBLIC
    ${LIBTAYLOR_INCLUDE_DIR}
  )

set_target_properties(pytaylor
  PROPERTIES
    PREFIX "${PYTHON_MODULE_PREFIX}"
    SUFFIX "${PYTHON_MODULE_EXTENSION}"
  )

file(
  COPY
    __init__.py
    version.py
  DESTINATION
    ${CMAKE_CURRENT_BINARY_DIR}
  )

install(
  TARGETS
    pytaylor
  LIBRARY DESTINATION .
  )
