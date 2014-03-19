# collect header files
SET (HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/AbstractDifferentialProblem.hpp 
	${CMAKE_CURRENT_SOURCE_DIR}/LinearDifferentialProblem.hpp
	${CMAKE_CURRENT_SOURCE_DIR}/NonLinearDifferentialProblem.hpp)

#  install target
#INSTALL (FILES ${HEADERS} DESTINATION "${HEADER_INSTALL_DIR}")
INSTALL_HEADERS_WITH_DIRECTORY ("${HEADERS}")
