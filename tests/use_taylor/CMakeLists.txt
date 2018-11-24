cmake_minimum_required(VERSION 3.5 FATAL_ERROR)

project(UseTaylor LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include(GNUInstallDirs)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_BINDIR})

find_package(taylor CONFIG REQUIRED)
get_property(_loc TARGET taylor::taylor PROPERTY INTERFACE_INCLUDE_DIRECTORIES)
message(STATUS "Found taylor: ${_loc} (found version ${taylor_VERSION})")

add_executable(use_taylor use_taylor.cpp)
target_link_libraries(use_taylor
  PUBLIC
  taylor::taylor
  )
