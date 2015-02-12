# collect header files
SET (HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/Expression.h
             ${CMAKE_CURRENT_SOURCE_DIR}/VariableExpression.h
             ${CMAKE_CURRENT_SOURCE_DIR}/DefaultEvaluator.h
             ${CMAKE_CURRENT_SOURCE_DIR}/TimeDependentExpression.h
             ${CMAKE_CURRENT_SOURCE_DIR}/expressions.h)

# collect source files
SET (SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/Expression.cpp
             ${CMAKE_CURRENT_SOURCE_DIR}/VariableExpression.cpp
             ${CMAKE_CURRENT_SOURCE_DIR}/DefaultEvaluator.cpp
             ${CMAKE_CURRENT_SOURCE_DIR}/TimeDependentExpression.cpp)

# add files to shared library dependancies
SET (CONTROL_PROBLEM_LIBRARY_SOURCES ${CONTROL_PROBLEM_LIBRARY_SOURCES} ${SOURCES} PARENT_SCOPE)
SET (CONTROL_PROBLEM_LIBRARY_HEADERS ${CONTROL_PROBLEM_LIBRARY_HEADERS} ${HEADERS} PARENT_SCOPE)

#  install headers
INSTALL_HEADERS_WITH_DIRECTORY ("${HEADERS}")
