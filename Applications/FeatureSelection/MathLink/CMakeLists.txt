# This is a SoPlex Wrapper directory for the project
# It should be built if you're using SoPlex

# Use our custom MPREP macro to custom-compile the .tm file
MMA_MPREP(
	${FSELECT_SOURCE_DIR}/MathLink/mmtemplate.tm
	${FSELECT_BINARY_DIR}/MathLink/mmtemplate.cpp)

ADD_EXECUTABLE(fselectML mlmain.cpp ${FSELECT_BINARY_DIR}/MathLink/mmtemplate.cpp)

TARGET_LINK_LIBRARIES(fselectML
	FSSoPlexWrapper
	FeatureSelection
	SoPlex
	${MMA_LIBRARIES}
    ${MYLIBS_BASE_LIBS}
    ITKNumerics
)

LINK_DIRECTORIES(
    ${FSELECT_BINARY_DIR}/FeatureSelection
    ${FSELECT_BINARY_DIR}/SoPlexWrapper
)