PROJECT(2DCMRep)

INCLUDE_DIRECTORIES(${2DCMRep_SOURCE_DIR}/Code)

ADD_EXECUTABLE(WaveletTest 
  Code/WaveletTest.cxx 
  Code/Wavelet.cxx)

ADD_EXECUTABLE(WaveletRepTest 
  Code/WaveletRepTest.cxx 
  Code/WaveletRep.cxx 
  Code/Wavelet.cxx)

ADD_EXECUTABLE(TemplatedWaveletRepTest 
  Code/TemplatedWaveletRepTest.cxx 
  Code/TemplatedWaveletRep.cxx 
  Code/WaveletRep.cxx 
  Code/Wavelet.cxx)

ADD_EXECUTABLE(MedialAtomTest 
  Code/MedialAtomTest.cxx 
  Code/MedialAtom.cxx)

ADD_EXECUTABLE(SparseMatrixTest 
  Code/SparseMatrixTest.cxx 
  Code/SparseMatrix.cxx)

ADD_EXECUTABLE(EquationSystemTest 
  Code/EquationSystemTest.cxx 
  Code/EquationSystem.cxx 
  Code/WaveletRep.cxx 
  Code/Wavelet.cxx)

ADD_EXECUTABLE(NewtonMethodTest 
  Code/NewtonMethodTest.cxx 
  Code/NewtonMethod.cxx 
  Code/EquationSystem.cxx 
  Code/PardisoInterface.cxx 
  Code/SparseMatrix.cxx 
  Code/WaveletRep.cxx 
  Code/Wavelet.cxx)

TARGET_LINK_LIBRARIES(NewtonMethodTest 
  ${PARDISO_LIB} 
  ${LAPACK_OR_MKL_LIBS} 
  ${PAULY_CORE_LIBS})

ADD_EXECUTABLE(CMRep2DTest 
  Code/CMRep2DTest.cxx  
  Code/CMRep2D.cxx 
  Code/NewtonMethod.cxx 
  Code/EquationSystem.cxx 
  Code/PardisoInterface.cxx 
  Code/SparseMatrix.cxx 
  Code/MedialAtom.cxx 
  Code/TemplatedWaveletRep.cxx 
  Code/WaveletRep.cxx 
  Code/Wavelet.cxx) 

TARGET_LINK_LIBRARIES(CMRep2DTest 
  ITKCommon 
  ITKIO 
  ${PARDISO_LIB} 
  ${LAPACK_OR_MKL_LIBS} 
  ${PAULY_CORE_LIBS})

ADD_EXECUTABLE(OptimizationTermsTest 
  Code/OptimizationTermsTest.cxx 
  Code/OptimizationTerms.cxx  
  Code/CMRep2D.cxx 
  Code/NewtonMethod.cxx 
  Code/EquationSystem.cxx 
  Code/PardisoInterface.cxx 
  Code/SparseMatrix.cxx 
  Code/MedialAtom.cxx 
  Code/TemplatedWaveletRep.cxx 
  Code/WaveletRep.cxx 
  Code/Wavelet.cxx) 

TARGET_LINK_LIBRARIES(OptimizationTermsTest 
  ITKCommon 
  ITKIO 
  ${PARDISO_LIB} 
  ${LAPACK_OR_MKL_LIBS} 
  ${PAULY_CORE_LIBS})

ADD_EXECUTABLE(ImageCheck Code/ImageCheck.cxx) 

TARGET_LINK_LIBRARIES(ImageCheck ITKCommon ITKIO )

