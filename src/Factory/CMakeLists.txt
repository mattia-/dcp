# create libraries
SET (SOURCES GenericFactory.hpp LinearSolverFactory.hpp LinearSolverProxy.cpp  Proxy.hpp)
SET (HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/GenericFactory.hpp 
	${CMAKE_CURRENT_SOURCE_DIR}/LinearSolverFactory.hpp 
	${CMAKE_CURRENT_SOURCE_DIR}/Proxy.hpp)

#ADD_LIBRARY (linearsolverproxy ${SOURCES}) # static library
ADD_LIBRARY (linearsolverproxy SHARED ${SOURCES}) # dynamic library

# INSTALLATION
INSTALL (TARGETS linearsolverproxy DESTINATION "${DYNAMIC_INSTALL_DIR}")
#INSTALL (FILES ${HEADERS} DESTINATION "${HEADER_INSTALL_DIR}")
INSTALL_HEADERS_WITH_DIRECTORY ("${HEADERS}")

