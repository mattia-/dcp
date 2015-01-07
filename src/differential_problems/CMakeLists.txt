# collect header files
SET (HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/AbstractDifferentialProblem.h 
             ${CMAKE_CURRENT_SOURCE_DIR}/LinearDifferentialProblem.h
             ${CMAKE_CURRENT_SOURCE_DIR}/NonlinearDifferentialProblem.h
             ${CMAKE_CURRENT_SOURCE_DIR}/TimeDependentDifferentialProblem.h
             ${CMAKE_CURRENT_SOURCE_DIR}/CompositeDifferentialProblem.h
             ${CMAKE_CURRENT_SOURCE_DIR}/SubdomainType.h
             ${CMAKE_CURRENT_SOURCE_DIR}/differential_problems.h)

# collect source files
SET (SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/CompositeDifferentialProblem.cpp
             ${CMAKE_CURRENT_SOURCE_DIR}/TimeDependentDifferentialProblem.cpp
             ${CMAKE_CURRENT_SOURCE_DIR}/AbstractDifferentialProblem.cpp)

# add files to shared library dependancies
SET (CONTROL_PROBLEM_LIBRARY_SOURCES ${CONTROL_PROBLEM_LIBRARY_SOURCES} ${SOURCES} PARENT_SCOPE)
SET (CONTROL_PROBLEM_LIBRARY_HEADERS ${CONTROL_PROBLEM_LIBRARY_HEADERS} ${HEADERS} PARENT_SCOPE)

#  install headers
INSTALL_HEADERS_WITH_DIRECTORY ("${HEADERS}")
