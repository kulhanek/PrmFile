# ==============================================================================
# PRMFile CMake File
# ==============================================================================

CMAKE_MINIMUM_REQUIRED(VERSION 2.5 FATAL_ERROR)
#SET(CMAKE_VERBOSE_MAKEFILE ON)

# ==============================================================================
# conditional setup ------------------------------------------------------------
# ==============================================================================

SET(LIBS_STATIC     OFF CACHE BOOL "Should the static version of prmfile library be built?")
SET(LIBS_SHARED     ON CACHE BOOL "Should the dynamic version of prmfile library be built?")
SET(TRY_QT_LIB      ON CACHE BOOL "Should the qt lib be used?")
SET(FORTRAN_ENABLED ON CACHE BOOL "Is Fortran90 available?")

# ==============================================================================
# project setup ----------------------------------------------------------------
# ==============================================================================

# define project name ----------------------------
SET(PRMFILE_VERSION 1.0.x)

# what languages we use --------------------------
ENABLE_LANGUAGE(CXX)
ENABLE_LANGUAGE(Fortran OPTIONAL)

IF(NOT CMAKE_Fortran_COMPILER_WORKS)
    SET(FORTRAN_ENABLED OFF)
ENDIF(NOT CMAKE_Fortran_COMPILER_WORKS)

# ==============================================================================
# external dependencies --------------------------------------------------------
# ==============================================================================

IF(NOT DEFINED COMPOSITE_PROJECT)
    SET(DEVELOPMENT_ROOT "$ENV{DEVELOPMENT_ROOT}")
ENDIF(NOT DEFINED COMPOSITE_PROJECT)

IF(NOT DEVELOPMENT_ROOT)
    MESSAGE(FATAL_ERROR "DEVELOPMENT_ROOT variable is not set!")
ENDIF(NOT DEVELOPMENT_ROOT)

# ==============================================================================
# external libraries -----------------------------------------------------------
# ==============================================================================

SET(HIPOLY_ROOT ${DEVELOPMENT_ROOT}/projects/hipoly/1.0)
INCLUDE_DIRECTORIES(${HIPOLY_ROOT}/src/lib/hipoly SYSTEM)
INCLUDE_DIRECTORIES(${HIPOLY_ROOT}/src/lib/hipoly/compress SYSTEM)
INCLUDE_DIRECTORIES(${HIPOLY_ROOT}/src/lib/hipoly/error SYSTEM)
INCLUDE_DIRECTORIES(${HIPOLY_ROOT}/src/lib/hipoly/event SYSTEM)
INCLUDE_DIRECTORIES(${HIPOLY_ROOT}/src/lib/hipoly/file SYSTEM)
INCLUDE_DIRECTORIES(${HIPOLY_ROOT}/src/lib/hipoly/list SYSTEM)
INCLUDE_DIRECTORIES(${HIPOLY_ROOT}/src/lib/hipoly/option SYSTEM)
INCLUDE_DIRECTORIES(${HIPOLY_ROOT}/src/lib/hipoly/string SYSTEM)
INCLUDE_DIRECTORIES(${HIPOLY_ROOT}/src/lib/hipoly/stream SYSTEM)
INCLUDE_DIRECTORIES(${HIPOLY_ROOT}/src/lib/hipoly/thread SYSTEM)
INCLUDE_DIRECTORIES(${HIPOLY_ROOT}/src/lib/hipoly/time SYSTEM)
INCLUDE_DIRECTORIES(${HIPOLY_ROOT}/src/lib/hipoly/xml SYSTEM)
LINK_DIRECTORIES(${HIPOLY_ROOT}/lib)
SET(HIPOLY_LIB_NAME hipoly)

# policies ---------------------------------------
IF(COMMAND cmake_policy)
    cmake_policy(SET CMP0003 NEW)
ENDIF(COMMAND cmake_policy)

# ==============================================================================
# output directories -----------------------------------------------------------
# ==============================================================================

SET(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR}/lib)

# ==============================================================================
# compiler setups --------------------------------------------------------------
# ==============================================================================

IF(NOT DEFINED COMPOSITE_PROJECT)
    # C++ options
    SET(CMAKE_CXX_FLAGS_RELEASE "-O2")
    SET(CMAKE_CXX_FLAGS_DEBUG "-g -Wall -Werror -pedantic-errors -Wundef -Wno-long-long")

    IF(FCOMP_NAME AND FORTRAN_ENABLED)
        # get the name of fortran compiler
        get_filename_component(FCOMP_NAME ${CMAKE_Fortran_COMPILER} NAME_WE)

        # default compiler setups - gfortran
        SET(CMAKE_Fortran_FLAGS_RELEASE "-O0")
        SET(CMAKE_Fortran_FLAGS_DEBUG   "-g -O0")

        # ifort
        IF(FCOMP_NAME STREQUAL ifort)
            SET(CMAKE_Fortran_FLAGS_RELEASE "-O2")
            SET(CMAKE_Fortran_FLAGS_DEBUG   "-g -O0 -std95 -traceback -check all -warn all")
        ENDIF(FCOMP_NAME STREQUAL ifort)

        # open64
        IF(FCOMP_NAME STREQUAL openf95)
            SET(CMAKE_Fortran_FLAGS_RELEASE "-O3 -fno-second-underscore")
            SET(CMAKE_Fortran_FLAGS_DEBUG "-g -O0 -fno-second-underscore")
        ENDIF(FCOMP_NAME STREQUAL openf95)

        # fix the problem with dynamic libraries and xlf compiler on SP5
        IF(FCOMP_NAME STREQUAL xlf)
            SET(CMAKE_Fortran_FLAGS_RELEASE "-O2 -qsuppress=1501-510:cmpmsg")
            SET(CMAKE_Fortran_FLAGS_DEBUG   "-g")
            IF(LIBS_SHARED)
                MESSAGE(FATAL_ERROR "LIBS_SHARED is ON, which is not supported with xlf compiler!")
            ENDIF(LIBS_SHARED)
        ENDIF(FCOMP_NAME STREQUAL xlf)

    ENDIF(FCOMP_NAME AND FORTRAN_ENABLED)

    SET(CMAKE_BUILD_TYPE Release)
    #SET(CMAKE_BUILD_TYPE Debug)
ENDIF(NOT DEFINED COMPOSITE_PROJECT)

# ==============================================================================
# project subdirectories  ------------------------------------------------------
# ==============================================================================

ADD_SUBDIRECTORY(src)
