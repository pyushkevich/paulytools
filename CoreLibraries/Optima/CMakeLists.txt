SET(OPTIMA_SRCS
    BrentLinearMethod.cpp
    ConjugateGradientMethod.cpp
    EvolutionaryStrategy.cpp
    optima.cpp
    Powell.cpp
    problems.cpp
    SimplexMethod.cpp
    support.cpp
)

ADD_LIBRARY(Optima ${OPTIMA_SRCS})
