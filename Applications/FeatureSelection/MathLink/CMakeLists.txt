# This is a SoPlex Wrapper directory for the project
# It should be built if you're using SoPlex

CONFIGURE_FILE(${FSELECT_SOURCE_DIR}/MathLink/mmtemplate.cpp.in ${FSELECT_BINARY_DIR}/MathLink/mmtemplate.cpp COPYONLY)

ADD_EXECUTABLE(fselectML mlmain.cpp ${FSELECT_BINARY_DIR}/MathLink/mmtemplate.cpp)

TARGET_LINK_LIBRARIES(fselectML
	FSSoPlexWrapper
	FeatureSelection
	${SOPLEX_LIB}
        ${MYLIBS_BASE_LIBS}
)

# Find Mathematica
IF(WIN32)
    FIND_PATH(MATHEMATICA_PATH Mathematica.exe 
        "C:/Program Files/Wolfram Research/Mathematica/4.1" 
        "D:/Program Files/Wolfram Research/Mathematica/4.1"
        "C:/Program Files/Wolfram Research/Mathematica/4.2" 
        "C:/Program Files/Wolfram Research/Mathematica/4.2" 
        "D:/Program Files/Wolfram Research/Mathematica/5.0")
        
    IF(MATHEMATICA_PATH)
        SET(MM_LIB_PATH "${MATHEMATICA_PATH}/AddOns/MathLink/DeveloperKit/Windows/CompilerAdditions/mldev32/lib")
        SET(MM_INCLUDE_PATH "${MATHEMATICA_PATH}/AddOns/MathLink/DeveloperKit/Windows/CompilerAdditions/mldev32/include")
        SET(MM_BIN_PATH "${MATHEMATICA_PATH}/AddOns/MathLink/DeveloperKit/Windows/CompilerAdditions/mldev32/bin")
		
        ADD_CUSTOM_COMMAND(
                SOURCE ${FSELECT_SOURCE_DIR}/MathLink/mmtemplate.tm
                COMMAND "${MM_BIN_PATH}/mprep.exe"
                ARGS ${FSELECT_SOURCE_DIR}/MathLink/mmtemplate.tm "-o" ${FSELECT_BINARY_DIR}/MathLink/mmtemplate.cpp
                OUTPUTS ${FSELECT_BINARY_DIR}/MathLink/mmtemplate.cpp
                TARGET fselectML)

        LINK_DIRECTORIES(${MM_LIB_PATH})
        INCLUDE_DIRECTORIES(${MM_INCLUDE_PATH})
        TARGET_LINK_LIBRARIES(fselectML ml32i1m ml32i2m)
    ENDIF(MATHEMATICA_PATH)
ENDIF(WIN32)

LINK_DIRECTORIES(
    ${FSELECT_BINARY_DIR}/FeatureSelection
    ${FSELECT_BINARY_DIR}/SoPlexWrapper
)

