PROJECT(FSELECT)

INCLUDE_DIRECTORIES(${FSELECT_SOURCE_DIR}/Generic)
INCLUDE_DIRECTORIES(${FSELECT_SOURCE_DIR}/SoPlexWrapper)
INCLUDE_DIRECTORIES(${FSELECT_SOURCE_DIR}/Standalone)

SUBDIRS(Generic SoPlexWrapper Standalone)

# Check if user wants to compile mathematica code
OPTION(USE_MATHEMATICA "Enable if you want to build MathLink interface for Mathematica")
IF(USE_MATHEMATICA)
    INCLUDE_DIRECTORIES(${FSELECT_SOURCE_DIR}/MathLink)    
    SUBDIRS(MathLink)
ENDIF(USE_MATHEMATICA)


