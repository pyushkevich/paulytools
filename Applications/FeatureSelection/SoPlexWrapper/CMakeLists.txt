# This is a SoPlex Wrapper directory for the project
# It should be built if you're using SoPlex

LINK_DIRECTORIES(${FSELECT_BINARY_DIR}/FeatureSelection)

ADD_LIBRARY(FSSoPlexWrapper SoPlexWrapper.cpp)


