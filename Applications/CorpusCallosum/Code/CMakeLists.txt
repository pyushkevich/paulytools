PROJECT(2DCMRep)

FIND_PACKAGE(ITK) 
IF(ITK_FOUND) 
	INCLUDE(${ITK_USE_FILE}) 
ELSE(ITK_FOUND) 
	MESSAGE(FATAL_ERROR "ITK not found. Please set ITK_DIR.") 
ENDIF(ITK_FOUND) 

SET(PAULY_SOURCE_DIR /mnt/data2/pauly/dev/paulytools)
SET(PAULY_INCLUDE_DIRS
  ${PAULY_SOURCE_DIR}/CoreLibraries/
  ${PAULY_SOURCE_DIR}/CoreLibraries/BasicImaging
  ${PAULY_SOURCE_DIR}/CoreLibraries/BSpline
  ${PAULY_SOURCE_DIR}/CoreLibraries/FLTKViewers
  ${PAULY_SOURCE_DIR}/CoreLibraries/ITKExtensions/WindowedSincFilter/
  ${PAULY_SOURCE_DIR}/CoreLibraries/Matrix
  ${PAULY_SOURCE_DIR}/CoreLibraries/MReps2d
  ${PAULY_SOURCE_DIR}/CoreLibraries/MReps2d/GLTools
  ${PAULY_SOURCE_DIR}/CoreLibraries/Optima
  ${PAULY_SOURCE_DIR}/CoreLibraries/Registry/
  ${PAULY_SOURCE_DIR}/ExternalLibraries/FLTKExtras/flIDE-0.3/src
  ${PAULY_SOURCE_DIR}/ExternalLibraries/SoPlex/soplex-1.2.1/src
)

INCLUDE_DIRECTORIES(${PAULY_INCLUDE_DIRS})
ADD_EXECUTABLE(test test.cxx)
ADD_EXECUTABLE(BSplineCubUniTest BSplineCubUniTest.cxx SumBSplineCubUni.cxx BSplineCubUni.cxx)
ADD_EXECUTABLE(WaveletTest WaveletTest.cxx Wavelet.cxx)
ADD_EXECUTABLE(WaveletRepTest WaveletRepTest.cxx WaveletRep.cxx Wavelet.cxx)
ADD_EXECUTABLE(TemplatedWaveletRepTest TemplatedWaveletRepTest.cxx TemplatedWaveletRep.cxx WaveletRep.cxx Wavelet.cxx)
ADD_EXECUTABLE(MedialAtomTest MedialAtomTest.cxx MedialAtom.cxx)
ADD_EXECUTABLE(SparseMatrixTest SparseMatrixTest.cxx SparseMatrix.cxx)
ADD_EXECUTABLE(EquationSystemTest EquationSystemTest.cxx EquationSystem.cxx WaveletRep.cxx Wavelet.cxx BSplineCubUni.cxx)
ADD_EXECUTABLE(NewtonMethodTest NewtonMethodTest.cxx NewtonMethod.cxx EquationSystem.cxx PardisoInterface.cxx SparseMatrix.cxx WaveletRep.cxx Wavelet.cxx BSplineCubUni.cxx)
TARGET_LINK_LIBRARIES(NewtonMethodTest /mnt/data2/pauly/opt/pardiso/lib/libpardiso_GNU_IA32.so /usr/lib/liblapack.so /usr/lib/libblas.so /usr/lib/libg2c.so.0)

ADD_EXECUTABLE(CMRep2DTest CMRep2DTest.cxx  CMRep2D.cxx NewtonMethod.cxx EquationSystem.cxx PardisoInterface.cxx SparseMatrix.cxx MedialAtom.cxx TemplatedWaveletRep.cxx WaveletRep.cxx Wavelet.cxx BSplineCubUni) 
TARGET_LINK_LIBRARIES(CMRep2DTest ITKCommon ITKIO  /mnt/data2/pauly/opt/pardiso/lib/libpardiso_GNU_IA32.so /usr/lib/liblapack.so /usr/lib/libblas.so /usr/lib/libg2c.so.0 /mnt/data2/pauly/dev/gccdbg/CoreLibraries/Matrix/libMatrix.a /mnt/data2/pauly/dev/gccdbg/CoreLibraries/Optima/libOptima.a)

ADD_EXECUTABLE(OptimizationTermsTest OptimizationTermsTest.cxx OptimizationTerms.cxx  CMRep2D.cxx NewtonMethod.cxx EquationSystem.cxx PardisoInterface.cxx SparseMatrix.cxx MedialAtom.cxx TemplatedWaveletRep.cxx WaveletRep.cxx Wavelet.cxx BSplineCubUni.cxx) 
TARGET_LINK_LIBRARIES(OptimizationTermsTest ITKCommon ITKIO  /mnt/data2/pauly/opt/pardiso/lib/libpardiso_GNU_IA32.so /usr/lib/liblapack.so /usr/lib/libblas.so /usr/lib/libg2c.so.0 /mnt/data2/pauly/dev/gccdbg/CoreLibraries/Matrix/libMatrix.a /mnt/data2/pauly/dev/gccdbg/CoreLibraries/Optima/libOptima.a)

ADD_EXECUTABLE(Optimization Optimization.cxx OptimizationTermsBSpline.cxx  CMRep2D.cxx NewtonMethod.cxx EquationSystem.cxx PardisoInterface.cxx SparseMatrix.cxx MedialAtom.cxx SumBSplineCubUni.cxx  BSplineCubUni.cxx) 
TARGET_LINK_LIBRARIES(Optimization ITKCommon ITKIO  /mnt/data2/pauly/opt/pardiso/lib/libpardiso_GNU_IA32.so /usr/lib/liblapack.so /usr/lib/libblas.so /usr/lib/libg2c.so.0 /mnt/data2/pauly/dev/gccdbg/CoreLibraries/Matrix/libMatrix.a /mnt/data2/pauly/dev/gccdbg/CoreLibraries/Optima/libOptima.a)

ADD_EXECUTABLE(ImageCheck ImageCheck.cxx) 
TARGET_LINK_LIBRARIES(ImageCheck ITKCommon ITKIO  /mnt/data2/pauly/opt/pardiso/lib/libpardiso_GNU_IA32.so /usr/lib/liblapack.so /usr/lib/libblas.so /usr/lib/libg2c.so.0 /mnt/data2/pauly/dev/gccdbg/CoreLibraries/Matrix/libMatrix.a /mnt/data2/pauly/dev/gccdbg/CoreLibraries/Optima/libOptima.a)

ADD_EXECUTABLE(DsmOptimization DsmOptimization.cxx OptimizationTermsBSpline.cxx  CMRep2D.cxx NewtonMethod.cxx EquationSystem.cxx PardisoInterface.cxx SparseMatrix.cxx MedialAtom.cxx SumBSplineCubUni.cxx  BSplineCubUni.cxx) 
TARGET_LINK_LIBRARIES(DsmOptimization ITKCommon ITKIO  /mnt/data2/pauly/opt/pardiso/lib/libpardiso_GNU_IA32.so /usr/lib/liblapack.so /usr/lib/libblas.so /usr/lib/libg2c.so.0 /mnt/data2/pauly/dev/gccdbg/CoreLibraries/Matrix/libMatrix.a /mnt/data2/pauly/dev/gccdbg/CoreLibraries/Optima/libOptima.a /mnt/data1/huiz/research/eclipse/workspace/Mala_3D/lib/libMalaNumerics.a)