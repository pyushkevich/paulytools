PROJECT(CMRepTool2d)

# For FLTK to link on Windows, add comctl32
IF(WIN32)
	LINK_LIBRARIES(comctl32 wsock32)
ENDIF(WIN32)

# Set up the sources
SET(SOURCES
        core.cpp
        dsltool.cpp
        ispace.cpp
        likehood.cpp
        misc.cpp
        readpgm.cpp
        reg.cpp
        scripting.cpp
        SplineMrep.cpp
        tools.cpp
        uibind.cpp
        undo.cpp
        windows.cpp)

# Set up the headers
SET(HEADERS
        array2d.h
        dispwin.h
        dsltool.h
        dslwin.h
        ispace.h
        likehood.h
        misc.h
        readpgm.h
        reg.h
        scripting.h
        SplineMrep.h
        tools.h
        uibind.h
        undo.h)

# Create an executable
ADD_EXECUTABLE(CMrepsTool2d ${SOURCES} ${HEADERS})

# Wrap the UI code
FLTK_WRAP_UI(SNAP ui.fl)

# Link with required libraries
LINK_LIBRARIES(
		mreps2d creps2d mrepsgl
    ${FFTW3_LIBRARY}
		${FLTKEXTRAS_LIB}
		${CIMAGE_LIB}
		${MYLIBS_BASE_LIBS} 
		${FLTK_LIBRARY} 
		${GLU_LIBRARY} 
		${OPENGL_LIBRARY})
