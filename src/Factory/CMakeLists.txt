# create libraries
SET (SOURCES GenericFactory.hpp LinearSolverFactory.hpp LinearSolverProxy.cpp  Proxy.hpp)
SET (HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/GenericFactory.hpp LinearSolverFactory.hpp Proxy.hpp)

#ADD_LIBRARY (linearsolverproxy ${SOURCES}) # static library
ADD_LIBRARY (linearsolverproxy SHARED ${SOURCES}) # dynamic library

# INSTALLATION
#  install target
INSTALL (TARGETS linearsolverproxy DESTINATION "${DYNAMIC_INSTALL_DIR}")
INSTALL (FILES ${HEADERS} DESTINATION "${HEADER_INSTALL_DIR}")

