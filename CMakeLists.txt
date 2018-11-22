cmake_minimum_required(VERSION 2.8)
project(libtaylor CXX)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
include_directories(${PROJECT_SOURCE_DIR})
add_executable(unittest_taylor ${PROJECT_SOURCE_DIR}/unittest_taylor.cpp)

enable_testing()
add_test(Taylor ${CMAKE_BINARY_DIR}/bin/unittest_taylor)

install(FILES micromath.hpp polymul.hpp taylor.hpp taylor_math.hpp
              ctaylor.hpp ctaylor_math.hpp tmath.hpp
	      DESTINATION include
	      )
