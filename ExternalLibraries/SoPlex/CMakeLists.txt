PROJECT(SOPLEX)

# The directory where the Soplex was extracted to 
SET(SOURCE_SUBDIR ${SOPLEX_SOURCE_DIR}/soplex-1.2.1/src)

SET(SOPLEX_LIB_SRCS 
    ${SOURCE_SUBDIR}/changesoplex.cpp
    ${SOURCE_SUBDIR}/didxset.cpp
    ${SOURCE_SUBDIR}/dsvector.cpp
    ${SOURCE_SUBDIR}/dvector.cpp
    ${SOURCE_SUBDIR}/enter.cpp
    ${SOURCE_SUBDIR}/factor.cpp
    ${SOURCE_SUBDIR}/forest.cpp
    ${SOURCE_SUBDIR}/idxset.cpp
    ${SOURCE_SUBDIR}/leave.cpp
    ${SOURCE_SUBDIR}/lpcolset.cpp
    ${SOURCE_SUBDIR}/lprow.cpp
    ${SOURCE_SUBDIR}/lprowset.cpp
    ${SOURCE_SUBDIR}/message.cpp
    ${SOURCE_SUBDIR}/mpsinput.cpp
    ${SOURCE_SUBDIR}/nameset.cpp
    ${SOURCE_SUBDIR}/slufactor.cpp
    ${SOURCE_SUBDIR}/solve.cpp
    ${SOURCE_SUBDIR}/soplex.cpp
    ${SOURCE_SUBDIR}/spxaggregatesm.cpp
    ${SOURCE_SUBDIR}/spxbasis.cpp
    ${SOURCE_SUBDIR}/spxbounds.cpp
    ${SOURCE_SUBDIR}/spxchangebasis.cpp
    ${SOURCE_SUBDIR}/spxdefaultpr.cpp
    ${SOURCE_SUBDIR}/spxdefaultrt.cpp
    ${SOURCE_SUBDIR}/spxdefines.cpp
    ${SOURCE_SUBDIR}/spxdesc.cpp
    ${SOURCE_SUBDIR}/spxdevexpr.cpp
    ${SOURCE_SUBDIR}/spxequilisc.cpp
    ${SOURCE_SUBDIR}/spxfastrt.cpp
    ${SOURCE_SUBDIR}/spxgeneralsm.cpp
    ${SOURCE_SUBDIR}/spxharrisrt.cpp
    ${SOURCE_SUBDIR}/spxhybridpr.cpp
    ${SOURCE_SUBDIR}/spxid.cpp
    ${SOURCE_SUBDIR}/spxio.cpp
    ${SOURCE_SUBDIR}/spxlp.cpp
    ${SOURCE_SUBDIR}/spxlpfread.cpp
    ${SOURCE_SUBDIR}/spxmpsread.cpp
    ${SOURCE_SUBDIR}/spxmpswrite.cpp
    ${SOURCE_SUBDIR}/spxparmultpr.cpp
    ${SOURCE_SUBDIR}/spxquality.cpp
    ${SOURCE_SUBDIR}/spxredundantsm.cpp
    ${SOURCE_SUBDIR}/spxrem1sm.cpp
    ${SOURCE_SUBDIR}/spxscaler.cpp
    ${SOURCE_SUBDIR}/spxshift.cpp
    ${SOURCE_SUBDIR}/spxsolve.cpp
    ${SOURCE_SUBDIR}/spxsolver.cpp
    ${SOURCE_SUBDIR}/spxstarter.cpp
    ${SOURCE_SUBDIR}/spxsteeppr.cpp
    ${SOURCE_SUBDIR}/spxsumst.cpp
    ${SOURCE_SUBDIR}/spxvecs.cpp
    ${SOURCE_SUBDIR}/spxvectorst.cpp
    ${SOURCE_SUBDIR}/spxweightpr.cpp
    ${SOURCE_SUBDIR}/spxweightst.cpp
    ${SOURCE_SUBDIR}/ssvector.cpp
    ${SOURCE_SUBDIR}/svector.cpp
    ${SOURCE_SUBDIR}/svset.cpp
    ${SOURCE_SUBDIR}/timer.cpp
    ${SOURCE_SUBDIR}/unitvector.cpp
    ${SOURCE_SUBDIR}/update.cpp
    ${SOURCE_SUBDIR}/updatevector.cpp
    ${SOURCE_SUBDIR}/vector.cpp
    ${SOURCE_SUBDIR}/vsolve.cpp
)

# On windows this generates some pragma warning!
IF(WIN32)
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /wd4274")
ENDIF(WIN32)

# Add long double option
OPTION(USE_LONG_DOUBLE "Compile with long double support")
IF(USE_LONG_DOUBLE)
	ADD_DEFINITIONS(-DWITH_LONG_DOUBLE)
ENDIF(USE_LONG_DOUBLE)

# Create the SoPlex library
ADD_LIBRARY(SoPlex ${SOPLEX_LIB_SRCS})

# Create the SoPlex example executable
ADD_EXECUTABLE(example ${SOURCE_SUBDIR}/example.cpp)
TARGET_LINK_LIBRARIES(example SoPlex)

# Set the include directory
SET(SOPLEX_INCLUDE_DIR ${SOURCE_SUBDIR})

# Set the library directory
SET(SOPLEX_LIBRARY_DIR ${SOPLEX_BINARY_DIR})

# Pass on the library and include paths to config file in binary directory
CONFIGURE_FILE(
	"${SOPLEX_SOURCE_DIR}/Config.cmake.in"
	"${SOPLEX_BINARY_DIR}/Config.cmake"
	IMMEDIATE @ONLY
)

