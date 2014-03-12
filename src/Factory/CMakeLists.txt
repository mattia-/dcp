# create libraries
SET (SOURCES GenericFactory.hpp LinearSolverFactory.hpp LinearSolverProxy.cpp  Proxy.hpp)
SET (HEADERS GenericFactory.hpp LinearSolverFactory.hpp Proxy.hpp)

#ADD_LIBRARY (linearsolverproxy ${SOURCES}) # static library
ADD_LIBRARY (linearsolverproxy SHARED ${SOURCES}) # dynamic library

# INSTALLATION
#  install target
INSTALL (TARGETS linearsolverproxy DESTINATION "${INSTALL_DIR}/lib")
INSTALL (FILES ${HEADERS} DESTINATION "${INSTALL_DIR}/include")
