# Mathematica CMake Search Module
IF(WIN32)
  FIND_PATH(MMA_ROOT_DIR Mathematica.exe 
    "C:/Program Files/Wolfram Research/Mathematica/4.1" 
    "D:/Program Files/Wolfram Research/Mathematica/4.1"
    "C:/Program Files/Wolfram Research/Mathematica/4.2" 
    "D:/Program Files/Wolfram Research/Mathematica/4.2"
    "C:/Program Files/Wolfram Research/Mathematica/5.0" 
    "D:/Program Files/Wolfram Research/Mathematica/5.0"
    )

  IF(MMA_ROOT_DIR)

    # Shorthand for the Mathlink directory
    SET(MMA_MATHLINK_DIR 
      "${MMA_ROOT_DIR}/AddOns/MathLink/DeveloperKit/Windows/CompilerAdditions/mldev32")

    # Set up the include directory
    FIND_PATH(MMA_INCLUDE_DIR mathlink.h "${MMA_MATHLINK_DIR}/include")

    # Set up the libraries
    FIND_LIBRARY(MMA_LIBRARY1 NAMES ml32i1m PATHS ${MMA_MATHLINK_DIR}/lib)
    FIND_LIBRARY(MMA_LIBRARY2 NAMES ml32i2m PATHS ${MMA_MATHLINK_DIR}/lib)
    SET(MMA_LIBRARIES ${MMA_LIBRARY1} ${MMA_LIBRARY2})

    # Set up the m-prep tool
    FIND_PROGRAM(MMA_MPREP_TOOL mprep 
      ${MMA_MATHLINK_DIR}/bin
      ${MMA_MATHLINK_DIR}/bin/*)

    IF(MMA_MPREP_TOOL)
      MACRO(MMA_MPREP ARG_TM ARG_CPP)
        ADD_CUSTOM_COMMAND(
          SOURCE ${ARG_TM}
          COMMAND ${MMA_MPREP_TOOL}
          ARGS ${ARG_TM} "-o" ${ARG_CPP}
          OUTPUT ${ARG_CPP})
      ENDMACRO(MMA_MPREP)
    ELSE(MMA_MPREP_TOOL)
      MESSAGE("Gotta find that tool [${MMA_MPREP_TOOL}]!")
    ENDIF(MMA_MPREP_TOOL)

  ENDIF(MMA_ROOT_DIR)

ENDIF(WIN32)

# Set the flag that MMA has been found
IF(MMA_LIBRARY1)
  IF(MMA_LIBRARY2)
    IF(MMA_MPREP_TOOL)
      SET (MMA_FOUND 1 CACHE INTERNAL "Mathematica and Mathlink are fully available!")
    ENDIF(MMA_MPREP_TOOL)
  ENDIF(MMA_LIBRARY2)
ENDIF(MMA_LIBRARY1)
