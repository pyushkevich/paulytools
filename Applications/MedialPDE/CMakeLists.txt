PROJECT(MedialPDE)

# Include the various directories
INCLUDE_DIRECTORIES(
  ${MedialPDE_SOURCE_DIR}/Code
  ${PAULY_SOURCE_DIR}/Applications/CMRepTool3d/Code)

# Sources for the PDE executable
SET(COMMON_SRCS
  Code/BasisFunctions2D.cxx
  Code/BranchingSubdivisionSurface.cxx
  Code/BruteForceSubdivisionMedialModel.cxx
  Code/CoefficientMapping.cxx
  Code/CartesianMedialModel.cxx
  Code/GeometryDescriptor.cxx
  Code/ITKImageWrapper.cxx
  Code/MedialAtom.cxx
  Code/MedialAtomGrid.cxx
  Code/MeshMedialPDESolver.cxx
  Code/MedialModelIO.cxx
  Code/MedialPDEMasks.cxx
  Code/MedialPDESites.cxx
  Code/MedialPDERenderer.cxx
  Code/OptimizationTerms.cxx
  Code/OptimizationParameters.cxx
  Code/PardisoInterface.cxx
  Code/PDESubdivisionMedialModel.cxx
  Code/PrincipalComponents.cxx
  Code/Procrustes.cxx
  Code/Registry.cxx
  Code/ScriptImaging.cxx
  Code/ScriptInterface.cxx
  Code/SparseMatrix.cxx
  Code/SubdivisionMedialModel.cxx
  Code/SubdivisionSurface.cxx
  Code/TestSolver.cxx
  Code/VTKMedialMesh.cxx
  Code/VTKSubdivision.cxx
  ${PAULY_SOURCE_DIR}/Applications/CMRepTool3d/Code/fracviewer.cpp
  ${PAULY_SOURCE_DIR}/Applications/CMRepTool3d/Code/glengine.cpp)

# Create a directory for SWIG output
FILE(MAKE_DIRECTORY ${MedialPDE_BINARY_DIR}/Code)

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
  vtkImaging
  vtkGraphics
  ${LAPACK_OR_MKL_LIBS}
  ${GSL_LIBRARY} ${GSL_CBLAS_LIBRARY}
  toms611
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

ADD_EXECUTABLE(testsub Code/TestSubdivision.cxx)
TARGET_LINK_LIBRARIES(testsub ${PDE_LIBS})
TARGET_LINK_LIBRARIES(testsub medialpde)

ADD_EXECUTABLE(cmr_fitbin Code/ProgramFitBinary.cxx ${COMMON_SRCS})
TARGET_LINK_LIBRARIES(cmr_fitbin ${PDE_LIBS})

ADD_EXECUTABLE(cmr_fit Code/Fit.cxx ${COMMON_SRCS})
TARGET_LINK_LIBRARIES(cmr_fit ${PDE_LIBS})

ADD_EXECUTABLE(cmr_sample_img Code/SampleImage.cxx ${COMMON_SRCS})
TARGET_LINK_LIBRARIES(cmr_sample_img ${PDE_LIBS})

# Create a statically linked profiling executable
### ADD_EXECUTABLE(mpdeprof Code/ProfilerDriver.cxx ${COMMON_SRCS})
### TARGET_LINK_LIBRARIES(mpdeprof ${PDE_LIBS} ${GLUT_LIBRARY})

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
# SET(MPDE_TEST_OBJMESH ${MPDE_TESTDATA_DIR}/triangle_mesh.obj)
SET(MPDE_TEST_OBJMESH ${MPDE_TESTDATA_DIR}/gr1a_gimpel_CN_snap_left_49pt.obj)
SET(MPDE_TEST_SUBMREP ${MPDE_TESTDATA_DIR}/gr1a_gimpel_CN_snap_left_49pt.cmrep)
SET(MPDE_TEST_BRUTE ${MPDE_TESTDATA_DIR}/cc_tempinit.cmrep)
SET(MPDE_TEST_BRUTE_IMG ${MPDE_TESTDATA_DIR}/cc_mask_blur)
SET(MPDE_TEST_BRANCH ${MPDE_TESTDATA_DIR}/tiny_branch.obj)

# Test various derivative computations
ADD_TEST(TestAffineTransform testpde AFFINE ${MPDE_TEST_CMREP})
ADD_TEST(TestVolumeComputation testpde VOLUME1 ${MPDE_TEST_CMREP})
ADD_TEST(TestDerivativesNoImage testpde DERIV1 ${MPDE_TEST_CMREP})
ADD_TEST(TestDerivativesWithImage testpde DERIV2 ${MPDE_TEST_CMREP})
ADD_TEST(TestBasisFunctionVariation testpde DERIV3 ${MPDE_TEST_CMREP})
ADD_TEST(TestDifferentialGeometry testpde DERIV4 ${MPDE_TEST_CMREP})
ADD_TEST(TestSparseMatrix testsub SPARSEMAT)
ADD_TEST(TestSubdivisionAffineTransform testpde AFFINE ${MPDE_TEST_SUBMREP})
ADD_TEST(TestSubdivisionSurface testsub SUBSURF1 ${MPDE_TEST_OBJMESH})
ADD_TEST(TestSubdivisionPDE testsub SUBSURF2 ${MPDE_TEST_OBJMESH})
ADD_TEST(TestSubSurfVolumeComputation testpde VOLUME1 ${MPDE_TEST_SUBMREP})
ADD_TEST(TestSubSurfDerivativesNoImage testpde DERIV1 ${MPDE_TEST_SUBMREP})
ADD_TEST(TestSubSurfDerivativesWithImage testpde DERIV2 ${MPDE_TEST_SUBMREP})
ADD_TEST(TestWedgeVolume testpde WEDGE)
ADD_TEST(TestModelSubdivision testsub MODELSUB ${MPDE_TEST_SUBMREP})
ADD_TEST(TestBranchingSubdivision testsub BRANCH1 ${MPDE_TEST_BRANCH})
ADD_TEST(TestBruteSubModelDerivativesNoImage testpde DERIV1 ${MPDE_TEST_BRUTE})
ADD_TEST(TestBruteSubModelDerivativesWithImage testpde DERIV2 ${MPDE_TEST_BRUTE} ${MPDE_TEST_BRUTE_IMG} vtk)



# Copy the script files
SET(SCRIPT_FILES
  cmrep_align_moments.py
  cmrep_brute2pde.py
  cmrep_build_images.py
  cmrep_cartesian_import.py
  cmrep_common.py
  cmrep_fitbinary.py
  cmrep_fit_discrete.py
  cmrep_getmesh.py
  cmrep_pca_compute.py
  cmrep_refine_fourier.py
  cmrep_render.py
  cmrep_sample_interior.py
  cmrep_subdivide.py
  )

FOREACH(FILE ${SCRIPT_FILES})
  CONFIGURE_FILE(
    ${MedialPDE_SOURCE_DIR}/Python/${FILE}
    ${MedialPDE_BINARY_DIR}/${FILE} )
ENDFOREACH(FILE ${SCRIPT_FILES})
