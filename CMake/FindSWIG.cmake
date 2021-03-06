# We must first find PYTHON
INCLUDE(${CMAKE_ROOT}/Modules/FindPythonLibs.cmake)

# Find the python binary
GET_FILENAME_COMPONENT(PYTHON_BIN_DIR ${PYHTON_INCLUDE_PATH}/../../bin ABSOLUTE)
MARK_AS_ADVANCED(PYTHON_BIN_DIR)
FIND_PROGRAM(PYTHON_BINARY python ${PYTHON_BIN_DIR})
  
# Search for Standard-SWIG (not ITK-based CMakeSwig)
FIND_PROGRAM(SWIG_PROGRAM swig)
IF(SWIG_PROGRAM)

  # Define the swig macro (python)
  MACRO(SWIG_WRAP_PYTHON_C ARG_INPUT ARG_OUTPUT_CPP)
    ADD_CUSTOM_COMMAND(
      SOURCE ${ARG_INPUT}
      COMMAND ${SWIG_PROGRAM}
      ARGS -python -o ${ARG_OUTPUT_CPP} ${ARG_INPUT}
      OUTPUT ${ARG_OUTPUT_CPP})
  ENDMACRO(SWIG_WRAP_PYTHON_C)

  # Define the swig macro (python, c++)
  MACRO(SWIG_WRAP_PYTHON_CXX ARG_INPUT ARG_OUTPUT_CPP)
    ADD_CUSTOM_COMMAND(
      SOURCE ${ARG_INPUT}
      COMMAND ${SWIG_PROGRAM}
      ARGS -python -c++ -o ${ARG_OUTPUT_CPP} ${ARG_INPUT}
      OUTPUT ${ARG_OUTPUT_CPP})
  ENDMACRO(SWIG_WRAP_PYTHON_CXX)

  # Define the swig macro (python, c++)
  MACRO(SWIG_WRAP_PYTHON_SHADOW ARG_INPUT ARG_OUTPUT_CPP ARG_OUTPUT_DIR ARG_MODULE)
    ADD_CUSTOM_COMMAND(
      SOURCE ${ARG_INPUT}
      COMMAND ${SWIG_PROGRAM}
      ARGS -python -shadow -c++ -o ${ARG_OUTPUT_CPP} -outdir ${ARG_OUTPUT_DIR} ${ARG_INPUT}
      OUTPUT ${ARG_OUTPUT_CPP})
  ENDMACRO(SWIG_WRAP_PYTHON_SHADOW)

ENDIF(SWIG_PROGRAM)
  
  

