cmake_minimum_required(VERSION 3.25...3.31)
project(
  KFParticle
  VERSION 26.1.21
  LANGUAGES CXX)

# prevent in-source builds #

if(PROJECT_SOURCE_DIR STREQUAL PROJECT_BINARY_DIR)
  message(FATAL_ERROR "You don't want to configure in the source directory!")
endif()

# C++ standard #

set(CMAKE_CXX_EXTENSIONS OFF) # ensure -std=c++xx instead of -std=g++xx
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_EXPORT_COMPILE_COMMANDS True)
set(CMAKE_CXX_SCAN_FOR_MODULES OFF) # remove -fmodules-ts flags

# includes and macros #

include(CheckIPOSupported)
include(CMakePackageConfigHelpers)
include(GNUInstallDirs) # GNU-like installation paths, e.g. lib/, include/, etc

# enable LTO #

check_ipo_supported(RESULT LTO_SUPPORTED)
if(LTO_SUPPORTED)
  message(STATUS "LTO/IPO is supported")
  set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)
else()
  message(STATUS "LTO/IPO is not supported")
endif()

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
target_compile_options(
  KFParticle PRIVATE -Wall -Wextra -Wshadow -Wnon-virtual-dtor -Wpedantic
                     -Wconversion)

# show available configs #

message("-- Available configurations: ${CMAKE_CONFIGURATION_TYPES}")

# apply config #

target_compile_options(
  KFParticle
  PRIVATE $<$<CONFIG:Release>:
          -O3
          -DNDEBUG
          -ffp-contract=fast
          -mtune=native
          -march=native>
          $<$<CONFIG:Debug>:
          -O0
          -g>)
target_compile_definitions(KFParticle PUBLIC $<$<CONFIG:Debug>:KF_DEBUG>)

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
