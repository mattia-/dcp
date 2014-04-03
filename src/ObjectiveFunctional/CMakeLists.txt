# collect header files
SET (HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/VariableExpression.hpp
             ${CMAKE_CURRENT_SOURCE_DIR}/AbstractObjectiveFunctional.hpp
             ${CMAKE_CURRENT_SOURCE_DIR}/ObjectiveFunctional.hpp)

# collect source files
SET (SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/VariableExpression.cpp
             ${CMAKE_CURRENT_SOURCE_DIR}/AbstractObjectiveFunctional.cpp)

# add files to shared library dependancies
SET (CONTROL_PROBLEM_LIBRARY_SOURCES ${CONTROL_PROBLEM_LIBRARY_SOURCES} ${SOURCES} PARENT_SCOPE)
SET (CONTROL_PROBLEM_LIBRARY_HEADERS ${CONTROL_PROBLEM_LIBRARY_HEADERS} ${HEADERS} PARENT_SCOPE)

#  install headers
#INSTALL (FILES ${HEADERS} DESTINATION "${HEADER_INSTALL_DIR}")
INSTALL_HEADERS_WITH_DIRECTORY ("${HEADERS}")
