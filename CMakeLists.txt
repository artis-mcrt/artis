cmake_minimum_required(VERSION 3.26)
project(artis)

set(default_build_type "Release")
# set(CMAKE_BUILD_TYPE Release)

set(CMAKE_BINARY_DIR ${CMAKE_SOURCE_DIR})
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})
set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR})

include_directories("${PROJECT_SOURCE_DIR}")

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

#set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)

set(SN3D_SOURCES sn3d.cc
decay.cc
gammapkt.cc
grid.cc
input.cc
kpkt.cc
ltepop.cc
macroatom.cc
md5.cc
nltepop.cc
nonthermal.cc
packet.cc
radfield.cc
ratecoeff.cc
rpkt.cc
spectrum_lightcurve.cc
stats.cc
thermalbalance.cc
update_grid.cc
update_packets.cc
vpkt.cc)

add_executable(sn3d ${SN3D_SOURCES})

set(EXSPEC_SOURCES exspec.cc
decay.cc
gammapkt.cc
grid.cc
input.cc
kpkt.cc
ltepop.cc
macroatom.cc
md5.cc
nltepop.cc
nonthermal.cc
packet.cc
radfield.cc
ratecoeff.cc
rpkt.cc
spectrum_lightcurve.cc
stats.cc
thermalbalance.cc
update_grid.cc
update_packets.cc
vpkt.cc)

add_executable(exspec ${EXSPEC_SOURCES})

find_package(MPI)
include_directories(SYSTEM ${MPI_INCLUDE_PATH})
target_link_libraries(sn3d MPI_CXX)

add_compile_options("-Wall" "-Wextra")

string(APPEND CMAKE_CXX_FLAGS_DEBUG "-g")
string(APPEND CMAKE_CXX_FLAGS_RELEASE "-O3")

#set(CMAKE_CXX_FLAGS "${CMAKE_C_FLAGS} -march=native")

add_compile_definitions(HAVE_INLINE)
add_compile_definitions(GSL_RANGE_CHECK_OFF)

# add_compile_options(-DHAVE_INLINE -DGSL_RANGE_CHECK_OFF)

target_link_libraries(sn3d gsl)
target_link_libraries(sn3d gslcblas)

target_link_libraries(exspec gsl)
target_link_libraries(exspec gslcblas)
