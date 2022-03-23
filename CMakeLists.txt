# Require a minimum version of cmake
cmake_minimum_required(VERSION 3.1)

project(rminimax C CXX)

message("Trying to build rminimax on " ${CMAKE_SYSTEM_NAME})

# Check that the we are using a C++11-conforming compiler
include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILRE_SUPPORTS_CXX0X)
if(COMPILER_SUPPORTS_CXX11)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
elseif(COMPILER_SUPPORTS_CXX0X)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
else()
    message(STATUS "The compiler ${CMAKE_CXX_COMPILER} has no C++11
    support. Please use a C++11-supporting compiler.")
endif()

# Check that the prerequisite libraries are installed and visible
message(STATUS "Starting dependency check...")
set (INCLUDE_DEPS gmp mpfr mpreal fplll qsopt_ex eigen3 flint2 arb)
set (LIB_DEPS gmp mpfr fplll qsopt_ex flint2 arb)

set (gmp_lib gmp)
set (mpfr_lib mpfr)
set (fplll_lib fplll)
set (qsopt_ex_lib qsopt_ex)
set (flint2_lib flint)
set (arb_lib arb)

set (gmp_header gmp.h)
set (mpfr_header mpfr.h)
set (mpreal_header mpreal.h)
set (fplll_header fplll.h)
set (qsopt_ex_header qsopt_ex/QSopt_ex.h)
set (eigen3_header eigen3/Eigen)
set (flint2_header flint/flint.h)
set (arb_header arb.h)

foreach (LIB ${LIB_DEPS})
    string (TOUPPER ${LIB} LIB_UPPER)
    find_library(${LIB_UPPER}_LIBRARY NAMES ${${LIB}_lib})
    if (NOT ${LIB_UPPER}_LIBRARY)
        message(FATAL_ERROR "\t ${LIB} library not found.")
    endif()
    add_library(${LIB} UNKNOWN IMPORTED)
    set_property(TARGET ${LIB} PROPERTY IMPORTED_LOCATION ${${LIB_UPPER}_LIBRARY})
    message("\t ${LIB} found in ${${LIB_UPPER}_LIBRARY}")

endforeach ()

foreach (LIB ${INCLUDE_DEPS})
    string(TOUPPER ${LIB} HEADER_PKG)
    set (HEADER ${${LIB}_header})
    find_path(${HEADER_PKG}_INCLUDE_DIR NAMES ${HEADER})
    if (NOT ${HEADER_PKG}_INCLUDE_DIR)
        message(FATAL_ERROR "\t ${HEADER} header not found.")
    endif()
    message("\t ${HEADER} found in ${${HEADER_PKG}_INCLUDE_DIR}")
    set (DEP_INCLUDE_DIRS ${DEP_INCLUDE_DIRS} ${${HEADER_PKG}_INCLUDE_DIR})
endforeach ()

message (STATUS "Finished dependency check.")

add_subdirectory(${PROJECT_SOURCE_DIR}/src)
add_definitions("-c")
set(COMMON_INCLUDES ${PROJECT_SOURCE_DIR}/rminimax)
include_directories(${DEP_INCLUDE_DIRS} ${COMMON_INCLUDES})

# constuct the main tool executable
add_executable(ratapprox ratapprox.cpp)
target_link_libraries(ratapprox ${LIB_DEPS} rminimax)
