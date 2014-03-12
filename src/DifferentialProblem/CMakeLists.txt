# collect header files
SET (HEADERS AbstractDifferentialProblem.hpp  LinearDifferentialProblem.hpp)

#  install target
INSTALL (FILES ${HEADERS} DESTINATION "${INSTALL_DIR}/include/")
