cmake_minimum_required(VERSION 3.0.0)
project(ptspostprocessing)
# c++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)
add_compile_options(-Wall -Wextra -pedantic -Werror)

# header and src files
include_directories(include)
file(GLOB SOURCES "src/*.cpp")

ADD_CUSTOM_TARGET(debug
  COMMAND ${CMAKE_COMMAND} -DCMAKE_BUILD_TYPE=Debug ${CMAKE_SOURCE_DIR}
  COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target all
  COMMENT "Switch CMAKE_BUILD_TYPE to Debug"
  )

ADD_CUSTOM_TARGET(release
  COMMAND ${CMAKE_COMMAND} -DCMAKE_BUILD_TYPE=Release ${CMAKE_SOURCE_DIR}
  COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target all
  COMMENT "Switch CMAKE_BUILD_TYPE to Release"
  )

# build shared library
add_library(LIBPROCESSING SHARED ${SOURCES})

# process executable
add_executable(process main.cpp)
target_link_libraries(process LIBPROCESSING)

# target
add_executable(checkcode test.cpp)
target_link_libraries(checkcode LIBPROCESSING)

# target
add_executable(dumppts dumppts.cpp)
target_link_libraries(dumppts LIBPROCESSING)

# target
add_executable(stabilitymode stabilitymode.cpp)
target_link_libraries(stabilitymode LIBPROCESSING)