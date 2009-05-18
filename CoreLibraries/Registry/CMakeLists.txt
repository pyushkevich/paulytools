SET(REGISTRY_SRCS
    Registry.cpp
    registry_old.cpp
    CommandLineArgumentParser.cxx
)

ADD_LIBRARY(Registry ${REGISTRY_SRCS})
