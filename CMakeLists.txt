cmake_minimum_required(VERSION 3.25...3.31)
project(KFParticle LANGUAGES CXX)

# prevent in-source builds #

if(PROJECT_SOURCE_DIR STREQUAL PROJECT_BINARY_DIR)
  message(FATAL_ERROR "You don't want to configure in the source directory!")
endif()

# c++ standard #

set(CMAKE_CXX_EXTENSIONS OFF) # ensure -std=c++xx instead of -std=g++xx
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_EXPORT_COMPILE_COMMANDS True)

# enable lto #

include(CheckIPOSupported)
check_ipo_supported(RESULT LTO_SUPPORTED)
if(LTO_SUPPORTED)
  message(STATUS "LTO/IPO is supported")
  set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)
else()
  message(STATUS "LTO/IPO is not supported")
endif()

# includes and macros #

include(CMakePackageConfigHelpers)
include(GNUInstallDirs)

# some definitions #

set(SOURCES KFParticle/KFParticle.cxx)
set(HEADERS KFParticle/KFParticle.hxx KFParticle/KFParticle_Math.hxx
            KFParticle/KFParticle_Const.hxx KFParticle/KFParticle_Utils.hxx)

# main target #

add_library(KFParticle SHARED ${SOURCES})
target_include_directories(
  KFParticle
  PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
         $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/KFParticle>
         $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)

# propagate c++23 to consumers
target_compile_features(KFParticle PUBLIC cxx_std_23)
target_compile_options(KFParticle PRIVATE -Wall -Wextra -Wshadow
                                          -Wnon-virtual-dtor)
target_compile_options(KFParticle PRIVATE -O3 -DNDEBUG -mtune=native
                                          -march=native -ffp-contract=fast)

option(KF_DEBUG "Print debug messages." OFF)
if(KF_DEBUG)
  target_compile_definitions(KFParticle PUBLIC KF_DEBUG)
endif(KF_DEBUG)

# installation of library and headers #

install(
  TARGETS KFParticle
  EXPORT ${PROJECT_NAME}Targets
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  INCLUDES
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

install(FILES ${HEADERS} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

# installation of package config #

configure_package_config_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/Config.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake
  INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}
  PATH_VARS CMAKE_INSTALL_INCLUDEDIR CMAKE_INSTALL_LIBDIR)

install(
  EXPORT ${PROJECT_NAME}Targets
  NAMESPACE ${PROJECT_NAME}::
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}
  EXPORT_LINK_INTERFACE_LIBRARIES)

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME})

# (optional) enable library from build tree #

export(
  EXPORT ${PROJECT_NAME}Targets
  NAMESPACE ${PROJECT_NAME}::
  FILE ${PROJECT_NAME}Targets.cmake)
