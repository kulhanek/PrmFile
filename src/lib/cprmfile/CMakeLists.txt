# ==============================================================================
# PRMFile CMake File
# ==============================================================================

# objects in library -----------------------------------------------------------
SET(CPRMFILE_SRC
        PrmFile.cpp
        PrmGroup.cpp
        PrmLine.cpp
        PrmSection.cpp
        PrmUtils.cpp
        )

# create static library --------------------------------------------------------
IF(LIBS_STATIC)
    ADD_LIBRARY(cprmfile_static STATIC ${CPRMFILE_SRC})
    SET_TARGET_PROPERTIES(cprmfile_static PROPERTIES
                            OUTPUT_NAME cprmfile
                            CLEAN_DIRECT_OUTPUT 1)
ENDIF(LIBS_STATIC)

# create shared library --------------------------------------------------------
IF(LIBS_SHARED)
    ADD_LIBRARY(cprmfile_shared SHARED ${CPRMFILE_SRC})
    TARGET_LINK_LIBRARIES(cprmfile_shared ${HIPOLY_LIB_NAME})
    SET_TARGET_PROPERTIES(cprmfile_shared PROPERTIES
                            OUTPUT_NAME cprmfile
                            CLEAN_DIRECT_OUTPUT 1
                            VERSION ${PRMFILE_VERSION})
    INSTALL(TARGETS cprmfile_shared
            DESTINATION lib)
ENDIF(LIBS_SHARED)
