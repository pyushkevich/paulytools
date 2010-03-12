PROJECT(Registration)


INCLUDE_DIRECTORIES(${PAULY_SOURCE_DIR}/Applications/MedialPDE/Code)
INCLUDE_DIRECTORIES(${PAULY_SOURCE_DIR}/CoreLibraries/ITKExtensions)

ADD_EXECUTABLE(evolreg EvolutionaryRegistration.cxx)

SET(ITK_LIBS
  BasicImaging
  ITKCommon
  ITKIO
  ITKVoxBoIO
  ITKBasicFilters
  ITKAlgorithms
  ITKNumerics
  ITKStatistics
)

SET(COMMON_LIBS
  ${ITK_LIBS}
)


TARGET_LINK_LIBRARIES(evolreg ${COMMON_LIBS})

INSTALL_TARGETS("/bin" evolreg)
