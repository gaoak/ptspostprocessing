cmake_minimum_required(VERSION 3.0.0)
project(ptspostprocessing)

# third party library
SET(FFTW_SEARCH_PATHS $ENV{LD_LIBRARY_PATH} $ENV{FFTW_HOME}/lib)
FIND_LIBRARY(FFTW_LIBRARY NAMES fftw3 fftw3f PATHS ${FFTW_SEARCH_PATHS})
FIND_PATH(FFTW_INCLUDE_DIR NAMES fftw3.h CACHE FILEPATH "FFTW include directory.")

link_directories(${FFTW_LIBRARY})
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
target_link_libraries(LIBPROCESSING  ${FFTW_LIBRARY})

# process executable
add_executable(process main.cpp)
target_link_libraries(process LIBPROCESSING)

# target
add_executable(checkcode test.cpp)
target_link_libraries(checkcode LIBPROCESSING)

# target
add_executable(dumppts dumppts.cpp)
target_link_libraries(dumppts LIBPROCESSING)

# target, baseflow + amplitude * mode
add_executable(stabilitymode stabilitymode.cpp)
target_link_libraries(stabilitymode LIBPROCESSING)

# target, baseflow + amplitude * mode
add_executable(combinefield combinefield.cpp)
target_link_libraries(combinefield LIBPROCESSING)

# target z-t plot
add_executable(hispts historypoints.cpp)
target_link_libraries(hispts LIBPROCESSING)

# target process vortex wake
add_executable(processwake VortexWake.cpp)
target_link_libraries(processwake LIBPROCESSING)

# target process LBM data
add_executable(processLBM processLBM.cpp)
target_link_libraries(processLBM LIBPROCESSING)

# target process Lagrangian points
add_executable(processLag processLag.cpp)
target_link_libraries(processLag LIBPROCESSING)

# target process channel flow
add_executable(processChannel ChannelFlow.cpp)
target_link_libraries(processChannel LIBPROCESSING)