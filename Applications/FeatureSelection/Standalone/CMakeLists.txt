# This is a SoPlex Wrapper directory for the project
# It should be built if you're using SoPlex

LINK_DIRECTORIES(
    ${FSELECT_BINARY_DIR}/FeatureSelection
    ${FSELECT_BINARY_DIR}/SoPlexWrapper
)

INCLUDE_DIRECTORIES(
    ${FSELECT_SOURCE_DIR}/SoPlexWrapper
)

ADD_EXECUTABLE(fselect main.cpp mathematicaio.cpp)

TARGET_LINK_LIBRARIES(fselect	
	FSSoPlexWrapper
	FeatureSelection
  SoPlex
  ITKNumerics
  ${PAULY_CORE_LIBS}
)
