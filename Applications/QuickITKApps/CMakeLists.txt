OPTION(USE_BOOST "Use BOOST graph library" ON) 

IF(USE_BOOST)
  INCLUDE_DIRECTORIES("/usr/local/include/boost-1_31")
ENDIF(USE_BOOST)

SUBDIRS(
  BrainAtlas
  AntiAlias
  FillSTLMesh
)
