OPTION(USE_BOOST "Use BOOST graph library" ON) 

IF(USE_BOOST)
  FIND_PATH( BOOST_INCLUDE_DIR "boost/config.hpp" "/usr/local/include/boost-1_31" "c:/boost")
  IF(BOOST_INCLUDE_DIR)
    INCLUDE_DIRECTORIES(${BOOST_INCLUDE_DIR})
  ENDIF(BOOST_INCLUDE_DIR)
ENDIF(USE_BOOST)

SUBDIRS(
  BrainAtlas
  AntiAlias
  FillSTLMesh
  Ovoids
)
