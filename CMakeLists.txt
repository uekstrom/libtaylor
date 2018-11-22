cmake_minimum_required(VERSION 3.5 FATAL_ERROR)
project(libtaylor
  LANGUAGES CXX
  VERSION 1.0.0
  )

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

include(GNUInstallDirs)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_BINDIR})

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
add_executable(unittest_taylor ${CMAKE_CURRENT_SOURCE_DIR}/unittest_taylor.cpp)

enable_testing()

add_test(
  NAME Taylor
  COMMAND $<TARGET_FILE:unittest_taylor>
  )

install(
  FILES
    ctaylor.hpp
    ctaylor_math.hpp
    micromath.hpp
    polymul.hpp
    taylor.hpp
    taylor_math.hpp
    tmath.hpp
  DESTINATION
    ${CMAKE_INSTALL_INCLUDEDIR}
  )
