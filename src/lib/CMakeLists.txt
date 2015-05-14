# ==============================================================================
# PRMFile CMake File
# ==============================================================================

# this is required for cprmfile library ----------------------------------------
INCLUDE_DIRECTORIES(cprmfile SYSTEM)

# include subdirectories -------------------------------------------------------
ADD_SUBDIRECTORY(cprmfile)

IF(FORTRAN_ENABLED)
    ADD_SUBDIRECTORY(fprmfile)
ENDIF(FORTRAN_ENABLED)
