cmake_minimum_required(VERSION 3.17 FATAL_ERROR)
project(SNAx VERSION 2.1 DESCRIPTION "Computational routines for axion and axion-like particle gamma-ray signatures from supernova SN1987A" LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 20)

include(GNUInstallDirs)
include(FetchContent)

if(CMAKE_SYSTEM_PROCESSOR MATCHES "x86_64")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
endif()
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")

# Check for Python 3 and PYBIND11 (if desired)
set(INTERNAL_PYBIND11_DIR "${PROJECT_SOURCE_DIR}/code/external/pybind11")
set(Python_FIND_VIRTUALENV FIRST)
find_package(Python3 COMPONENTS Interpreter Development)
find_package(pybind11 CONFIG QUIET)
if(PYTHON3_FOUND)
  if(pybind11_FOUND)
    message("-- Python 3 and associated pybind11 found!")
  elseif(EXISTS ${INTERNAL_PYBIND11_DIR})
    message("-- Python 3 and internal pybind11 directory found!")
    add_subdirectory(${INTERNAL_PYBIND11_DIR})
  else()
    FetchContent_Declare(pybind11
      GIT_REPOSITORY https://github.com/pybind/pybind11.git
      GIT_TAG v2.9.2
      GIT_PROGRESS ON
      GIT_SHALLOW ON
      GIT_SUBMODULES ""
      GIT_CONFIG advice.detachedHead=False
      SOURCE_DIR ${INTERNAL_PYBIND11_DIR}
      )
    message("-- Python 3 found but not pybind11. Attempting to clone pybind11 from Github...")
    FetchContent_Populate(pybind11)
    add_subdirectory(${INTERNAL_PYBIND11_DIR})
  endif()
else()
  message(FATAL_ERROR "-- Could not find the Python Interpreter! Please provide provide a hint to CMAKE, using '-D CMAKE_PYTHON_EXECUTABLE'.")
endif()

add_subdirectory(code/cpp/)
