# This directory generates a generic library that need an LP wrapper to operate
ADD_LIBRARY(
  FeatureSelection 
  FeatureSelection.cpp 
  LinearSeparation.cpp 
  LPWrapper.cpp 
  WindowSelection.cpp 
  CrossValidation.cpp 
  globals.cpp
)


