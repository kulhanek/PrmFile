# ==============================================================================
# PRMFile CMake File
# ==============================================================================

# objects in library -----------------------------------------------------------
SET(FPRMFILE_OBJ
        prmfile_core.f90
        prmfile_dat.f90
        prmfile.F90
        )


# create static library --------------------------------------------------------
IF(LIBS_STATIC)
    ADD_LIBRARY(fprmfile_static STATIC ${FPRMFILE_OBJ})
    SET_TARGET_PROPERTIES(fprmfile_static PROPERTIES
                            OUTPUT_NAME fprmfile
                            CLEAN_DIRECT_OUTPUT 1)
ENDIF(LIBS_STATIC)

# create shared library --------------------------------------------------------
IF(LIBS_SHARED)
    ADD_LIBRARY(fprmfile_shared SHARED ${FPRMFILE_OBJ})
    SET_TARGET_PROPERTIES(fprmfile_shared PROPERTIES
                            OUTPUT_NAME fprmfile
                            CLEAN_DIRECT_OUTPUT 1
                            VERSION ${PRMFILE_VERSION})
    INSTALL(TARGETS fprmfile_shared
            DESTINATION lib)
ENDIF(LIBS_SHARED)

