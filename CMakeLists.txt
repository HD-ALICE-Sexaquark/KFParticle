cmake_minimum_required(VERSION 3.25...3.31)
project(
  KFParticle
  VERSION 26.1.28
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

# find dependency #

find_package(Armadillo REQUIRED CONFIG)

# main target #

add_library(${PROJECT_NAME} SHARED ${SOURCES})

target_include_directories(
  ${PROJECT_NAME}
  PRIVATE ${ARMADILLO_INCLUDE_DIRS}
          $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
          $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/KFParticle>
          $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)

target_link_libraries(${PROJECT_NAME} PRIVATE ${ARMADILLO_LIBRARIES})

target_compile_options(
  ${PROJECT_NAME} PRIVATE -Wall -Wextra -Wshadow -Wnon-virtual-dtor -Wpedantic
                          -Wconversion)

# set KF_DEBUG #

target_compile_definitions(${PROJECT_NAME} PUBLIC $<$<CONFIG:Debug>:KF_DEBUG>)

# installation of library and headers #

install(
  TARGETS ${PROJECT_NAME}
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
