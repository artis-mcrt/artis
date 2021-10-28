cmake_minimum_required(VERSION 3.9)
project(artis)

set(WALLTIMEHOURS 24)

set(default_build_type "Release")
# set(CMAKE_BUILD_TYPE Release)

set(CMAKE_BINARY_DIR ${CMAKE_SOURCE_DIR})
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})
set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR})

include_directories("${PROJECT_SOURCE_DIR}")

set(CMAKE_CXX_STANDARD 17)

#check_ipo_supported()
set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)

set(SN3D_SOURCES sn3d.cc atomic.cc boundary.cc emissivities.cc gamma.cc globals.cc grey_emissivities.cc grid.cc input.cc kpkt.cc light_curve.cc ltepop.cc macroatom.cc nltepop.cc nonthermal.cc decay.cc packet_init.cc photo_electric.cc polarization.cc radfield.cc ratecoeff.cc rpkt.cc stats.cc thermalbalance.cc update_grid.cc update_packets.cc vectors.cc vpkt.cc md5.cc
add_executable(sn3d ${SN3D_SOURCES})

set(EXSPEC_SOURCES exspec.cc grid.cc globals.cc input.cc vectors.cc packet_init.cc update_grid.cc update_packets.cc gamma.cc boundary.cc macroatom.cc decay.cc rpkt.cc kpkt.cc photo_electric.cc emissivities.cc grey_emissivities.cc ltepop.cc atomic.cc ratecoeff.cc thermalbalance.cc light_curve.cc spectrum.cc polarization.cc nltepop.cc nonthermal.cc radfield.cc stats.cc vpkt.cc md5.cc)

add_executable(exspec ${EXSPEC_SOURCES})
add_executable(exgamma ${EXSPEC_SOURCES})

if(UNIX AND NOT APPLE)
  find_package(MPI)
  include_directories(SYSTEM ${MPI_INCLUDE_PATH})
  target_link_libraries(sn3d MPI_CXX)
endif()

add_compile_options("-Wall" "-Wextra")

string(APPEND CMAKE_CXX_FLAGS_DEBUG "-g")
string(APPEND CMAKE_CXX_FLAGS_RELEASE "-O3")

set(CMAKE_CXX_FLAGS "${CMAKE_C_FLAGS} -march=native")

add_compile_definitions(HAVE_INLINE)
add_compile_definitions(GSL_RANGE_CHECK_OFF)

# add_compile_options(-DHAVE_INLINE -DGSL_RANGE_CHECK_OFF)

target_compile_definitions(sn3d PRIVATE WALLTIMELIMITSECONDS=${WALLTIMEHOURS}*3600)
target_link_libraries(sn3d gsl)
target_link_libraries(sn3d gslcblas)

target_link_libraries(exspec gsl)
target_link_libraries(exspec gslcblas)

target_link_libraries(exgamma gsl)
target_link_libraries(exgamma gslcblas)
