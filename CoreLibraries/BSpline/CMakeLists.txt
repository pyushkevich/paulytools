IF(BUILD_TESTING)
  ADD_EXECUTABLE(BSplineCurveTest BSplineCurveTest.cxx  BSplineCurve.h  BSplineCurve.txx)
  LINK_LIBRARIES(ITKNumerics)
ENDIF(BUILD_TESTING)
