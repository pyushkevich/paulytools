OPTION(BUILD_APP_CMRepTool2d "Build the tool for 2D mrep/cmrep fitting" OFF)
OPTION(BUILD_APP_CMRepTool3d "Build the tool for 3D mrep/cmrep fitting" OFF)
OPTION(BUILD_APP_FeatureSelection "Build the feature selection toolkit" OFF)
OPTION(BUILD_APP_QuickITKApps "Build small ITK applications (ITK QuickApps)" OFF)
OPTION(BUILD_APP_GraphPartition "Build graph partitioning code" OFF)
OPTION(BUILD_APP_MedialPDE "Build the PDE CM-Reps code" OFF)
OPTION(BUILD_APP_BrainTracer "Build the BrainTracer Application" OFF)

VARIABLE_REQUIRES(BUILD_APP_CMRepTool2d DO_BUILD_CMRepTool2d 
  BUILD_EXTERNAL
  USE_OPENGL
  USE_FLTK
  ITK_FOUND
)

VARIABLE_REQUIRES(BUILD_APP_CMRepTool3d DO_BUILD_CMRepTool3d 
  USE_OPENGL
  USE_GLUT
  ITK_FOUND
)

VARIABLE_REQUIRES(BUILD_APP_MedialPDE DO_BUILD_MedialPDE 
  USE_OPENGL
  USE_SWIG 
  USE_GLUT
  VTK_FOUND
  ITK_FOUND
  LAPACK_OR_MKL_FOUND
  PARDISO_FOUND
) 

VARIABLE_REQUIRES(BUILD_APP_FeatureSelection DO_BUILD_FeatureSelection
  BUILD_EXTERNAL
  USE_ITK
)

VARIABLE_REQUIRES(BUILD_APP_QuickITKApps DO_BUILD_QuickITKApps
  USE_ITK
  USE_VTK
  USE_OPENGL
  USE_FLTK
)

IF(BUILD_APP_CMRepTool2d)
  SUBDIRS(CMRepTool2d)
ENDIF(BUILD_APP_CMRepTool2d)

IF(BUILD_APP_CMRepTool3d)
  SUBDIRS(CMRepTool3d)
ENDIF(BUILD_APP_CMRepTool3d)

IF(BUILD_APP_FeatureSelection)
  SUBDIRS(FeatureSelection)
ENDIF(BUILD_APP_FeatureSelection)

IF(BUILD_APP_QuickITKApps)
  SUBDIRS(QuickITKApps)
ENDIF(BUILD_APP_QuickITKApps)

IF(BUILD_APP_GraphPartition)
  SUBDIRS(GraphPartitioning)
ENDIF(BUILD_APP_GraphPartition)

IF(BUILD_APP_MedialPDE)
  SUBDIRS(MedialPDE)
ENDIF(BUILD_APP_MedialPDE)

IF(BUILD_APP_BrainTracer)
  SUBDIRS(BrainTracer)
ENDIF(BUILD_APP_BrainTracer)

