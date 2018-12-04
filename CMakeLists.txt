cmake_minimum_required(VERSION 3.9)
project(artis)

set(WALLTIMEHOURS 24)

set(default_build_type "Release")
# set(CMAKE_BUILD_TYPE Release)

set(CMAKE_BINARY_DIR ${CMAKE_SOURCE_DIR})
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})
set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR})

include_directories("${PROJECT_SOURCE_DIR}")

set(CMAKE_C_STANDARD 11)

#check_ipo_supported()
set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)

set(SN3D_SOURCES sn3d.c atomic.c boundary.c compton.c emissivities.c gamma.c grey_emissivities.c grid_init.c input.c kpkt.c ltepop.c macroatom.c move.c nltepop.c nonthermal.c packet_init.c photo_electric.c polarization.c radfield.c ratecoeff.c rpkt.c thermalbalance.c time_init.c update_grid.c update_packets.c vectors.c vpkt.c md5.c)
add_executable(sn3d ${SN3D_SOURCES})

set(EXSPEC_SOURCES exspec.c grid_init.c input.c vectors.c packet_init.c time_init.c update_grid.c update_packets.c gamma.c boundary.c move.c compton.c macroatom.c rpkt.c kpkt.c photo_electric.c emissivities.c grey_emissivities.c ltepop.c atomic.c ratecoeff.c thermalbalance.c light_curve.c spectrum.c polarization.c nltepop.c radfield.c nonthermal.c vpkt.c md5.c)
add_executable(exspec ${EXSPEC_SOURCES})

set(EXGAMMA_SOURCES exgamma.c grid_init.c input.c vectors.c packet_init.c time_init.c update_grid.c update_packets.c gamma.c boundary.c move.c compton.c macroatom.c rpkt.c kpkt.c photo_electric.c emissivities.c grey_emissivities.c ltepop.c atomic.c ratecoeff.c thermalbalance.c light_curve.c spectrum.c polarization.c nltepop.c radfield.c nonthermal.c vpkt.c md5.c)
add_executable(exgamma ${EXGAMMA_SOURCES})

if(UNIX AND NOT APPLE)
  find_package(MPI)
  include_directories(SYSTEM ${MPI_INCLUDE_PATH})
  target_link_libraries(sn3d MPI_C)
endif()

add_compile_options("-Wall" "-Wextra")

string(APPEND CMAKE_C_FLAGS_DEBUG "-g")
string(APPEND CMAKE_C_FLAGS_RELEASE "-O3")

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -march=native")

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
