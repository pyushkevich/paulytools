OPTION(BUILD_CMRepTool2d "Build the tool for 2D mrep/cmrep fitting" OFF)
OPTION(BUILD_CMRepTool3d "Build the tool for 3D mrep/cmrep fitting" OFF)
OPTION(BUILD_FeatureSelection "Build the feature selection toolkit" OFF)

VARIABLE_REQUIRES(BUILD_CMRepTool2d DO_BUILD_CMRepTool2d 
  BUILD_EXTERNAL
  USE_OPENGL
  USE_FLTK
  USE_ITK
)

VARIABLE_REQUIRES(BUILD_CMRepTool2d DO_BUILD_CMRepTool3d 
  BUILD_EXTERNAL
  USE_OPENGL
  USE_FLTK
  USE_ITK
)

VARIABLE_REQUIRES(BUILD_FeatureSelection DO_BUILD_FeatureSelection
  BUILD_EXTERNAL
  USE_ITK
)

IF(BUILD_CMRepTool2d)
  SUBDIRS(CMRepTool2d)
ENDIF(BUILD_CMRepTool2d)

IF(BUILD_CMRepTool3d)
  SUBDIRS(CMRepTool3d)
ENDIF(BUILD_CMRepTool3d)

IF(BUILD_FeatureSelection)
  SUBDIRS(FeatureSelection)
ENDIF(BUILD_FeatureSelection)


