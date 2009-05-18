SET(REGISTRY_SRCS
    Registry.cxx
    registry_old.cpp
    CommandLineArgumentParser.cxx
)

ADD_LIBRARY(Registry ${REGISTRY_SRCS})
