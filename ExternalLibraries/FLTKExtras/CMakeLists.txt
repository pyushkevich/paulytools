PROJECT(FLTKEXTRAS)

SET(SOURCES
        flIDE-0.3/src/Fl_Node.cxx
        flIDE-0.3/src/Fl_ToggleNode.cxx
        flIDE-0.3/src/Fl_ToggleTree.cxx
        flIDE-0.3/src/Fl_Tree.cxx
)

# Set up FLTK
INCLUDE (${CMAKE_ROOT}/Modules/FindFLTK.cmake)

INCLUDE_DIRECTORIES(
        ${FLTKEXTRAS_SOURCE_DIR}/flIDE-0.3/src
        ${FLTK_INCLUDE_PATH}
)

ADD_LIBRARY(fltkextras ${SOURCES})

# Pass on the library and include paths to config file in binary directory
SET(FLTKEXTRAS_INCLUDE_DIR ${FLTKEXTRAS_SOURCE_DIR}/flIDE-0.3/src)
SET(FLTKEXTRAS_LIBRARY_DIR ${FLTKEXTRAS_BINARY_DIR})
CONFIGURE_FILE(
	"${FLTKEXTRAS_SOURCE_DIR}/Config.cmake.in"
	"${FLTKEXTRAS_BINARY_DIR}/Config.cmake"
	IMMEDIATE @ONLY
)
