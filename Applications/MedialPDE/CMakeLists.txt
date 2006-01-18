PROJECT(MedialPDE)

# Include the various directories
INCLUDE_DIRECTORIES(
  ${MedialPDE_SOURCE_DIR}/Code 
  ${PAULY_SOURCE_DIR}/Applications/CMRepTool3d/Code)

# Sources for the PDE executable
SET(COMMON_SRCS 
  ${PAULY_SOURCE_DIR}/Applications/CMRepTool3d/Code/fracviewer.cpp
  ${PAULY_SOURCE_DIR}/Applications/CMRepTool3d/Code/glengine.cpp
  Code/BasisFunctions2D.cxx
  Code/CartesianMedialAtomGrid.cxx
  Code/CoefficientMask.cxx
  Code/GeometryDescriptor.cxx
  Code/ITKImageWrapper.cxx
  Code/MedialAtom.cxx
  Code/MedialAtomGrid.cxx
  Code/MedialPDERenderer.cxx
  Code/MedialPDEMasks.cxx 
  Code/MedialPDESites.cxx 
  Code/MedialPDESolver.cxx 
  Code/OptimizationTerms.cxx
  Code/PardisoInterface.cxx
  Code/PrincipalComponents.cxx
  Code/Procrustes.cxx
  Code/Registry.cxx
  Code/ScriptImaging.cxx
  Code/ScriptInterface.cxx
  Code/TestSolver.cxx
  Code/VTKMedialMesh.cxx)

# Apply SWIG to the .i file
SET(SWIG_SOURCE ${MedialPDE_SOURCE_DIR}/Code/medialpde.i)
SET(SWIG_DEPENDS ${SWIG_SOURCE}  ${MedialPDE_SOURCE_DIR}/Code/ScriptInterface.h)
SET(SWIG_TARGET ${MedialPDE_BINARY_DIR}/Code/ScriptInterface_wrap.cxx)

ADD_CUSTOM_COMMAND(
  SOURCE ${SWIG_SOURCE} 
  COMMAND ${SWIG_PROGRAM}
  ARGS -python -shadow -c++ -o ${SWIG_TARGET} -outdir ${MedialPDE_BINARY_DIR} ${SWIG_SOURCE}
  DEPENDS ${SWIG_DEPENDS} 
  OUTPUT ${SWIG_TARGET}
  COMMENT "Python wrappers with SWIG")

# Add python to the includes
INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_PATH})

# The sources for the script-wrapped library
SET(LIB_SRCS ${COMMON_SRCS}
  ${MedialPDE_BINARY_DIR}/Code/ScriptInterface_wrap.cxx)

# Set the library set
SET(PDE_LIBS
  ${PAULY_CORE_LIBS}
  ${PARDISO_LIB} 
  ITKNumerics
  ITKIO
  ITKCommon
  ITKBasicFilters
  ITKVoxBoIO
  ${GLUT_LIBRARY}
  vtkCommon
  vtkIO
  ${LAPACK_OR_MKL_LIBS}
  ${GSL_LIBRARY} ${GSL_CBLAS_LIBRARY} 
  )

# Create a common library
ADD_LIBRARY(medialpde SHARED ${COMMON_SRCS})
TARGET_LINK_LIBRARIES(medialpde ${PDE_LIBS})

# Create the Python module
ADD_LIBRARY("_medialpdemodule" SHARED ${SWIG_TARGET})
TARGET_LINK_LIBRARIES("_medialpdemodule" medialpde)
SET_TARGET_PROPERTIES("_medialpdemodule" PROPERTIES PREFIX "")

# Create the test executable
ADD_EXECUTABLE(testpde Code/ScriptTest.cxx)
TARGET_LINK_LIBRARIES(testpde ${PDE_LIBS})
TARGET_LINK_LIBRARIES(testpde medialpde)

# Specify the data root
FIND_PATH(MPDE_DATA_ROOT "t1")
IF(MPDE_DATA_ROOT)
  FIND_PATH(MPDE_INPUTDATA_ROOT "init/surf01.txt" ${MPDE_DATA_ROOT}/Input)
  FIND_PATH(MPDE_TEMP_ROOT "" ${MPDE_DATA_ROOT}/tmp "/tmp")
ENDIF(MPDE_DATA_ROOT)
  

# Generate tests
ENABLE_TESTING()
SET(MPDE_TESTDATA_DIR ${MedialPDE_SOURCE_DIR}/Testing/Input)
SET(MPDE_TEST_CMREP ${MPDE_TESTDATA_DIR}/st1001L.align.mpde)
SET(MPDE_TEST_PCA ${MPDE_TESTDATA_DIR}/shapemat.mat)

# Test various derivative computations
ADD_TEST(TestDerivativesNoImage testpde 
  DERIV1 ${MPDE_TEST_CMREP} ${MPDE_TEST_PCA})

ADD_TEST(TestDerivativesWithImage testpde DERIV2 ${MPDE_TEST_CMREP})
ADD_TEST(TestBasisFunctionVariation testpde DERIV3 ${MPDE_TEST_CMREP})
ADD_TEST(TestDifferentialGeometry testpde DERIV4 ${MPDE_TEST_CMREP})

# Copy the script files
SET(SCRIPT_FILES 
  blurimage.py
  example.py 
  example2.py 
  example3.py 
  example4.py
  example5.py
  example6.py
  common.py 
  common2.py 
  imageprep.py
  segment.py
  segment2.py
  segment2.py
  segment3.py
  pipeline.py
  pipeline2.py
  pca.py
  pcalite.py
  pca_segment.py
  sample.py
  sample_img.py
  sample_ref.py
  testdphi.py
  vis.py
  vis_sample.py
  )

FOREACH(FILE ${SCRIPT_FILES})  
  CONFIGURE_FILE(
    ${MedialPDE_SOURCE_DIR}/Python/${FILE}
    ${MedialPDE_BINARY_DIR}/${FILE} )
ENDFOREACH(FILE ${SCRIPT_FILES})  
