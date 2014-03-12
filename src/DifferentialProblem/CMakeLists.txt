# collect header files
SET (HEADERS AbstractDifferentialProblem.hpp  LinearDifferentialProblem.hpp)

#  install target
INSTALL (FILES ${HEADERS} DESTINATION "${HEADER_INSTALL_DIR}")
