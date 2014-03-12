# collect header files
SET (HEADERS AbstractDifferentialProblem.hpp  LinearDifferentialProblem.hpp)

#  install target
INSTALL (TARGETS ${HEADERS} "${INSTALL_DIR}/include/")
