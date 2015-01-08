# collect header files
SET (HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/AbstractProblem.h 
             ${CMAKE_CURRENT_SOURCE_DIR}/LinearProblem.h
             ${CMAKE_CURRENT_SOURCE_DIR}/NonlinearProblem.h
             ${CMAKE_CURRENT_SOURCE_DIR}/TimeDependentProblem.h
             ${CMAKE_CURRENT_SOURCE_DIR}/CompositeProblem.h
             ${CMAKE_CURRENT_SOURCE_DIR}/SubdomainType.h
             ${CMAKE_CURRENT_SOURCE_DIR}/differential_problems.h)

# collect source files
SET (SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/CompositeProblem.cpp
             ${CMAKE_CURRENT_SOURCE_DIR}/TimeDependentProblem.cpp
             ${CMAKE_CURRENT_SOURCE_DIR}/AbstractProblem.cpp)

# add files to shared library dependancies
SET (CONTROL_PROBLEM_LIBRARY_SOURCES ${CONTROL_PROBLEM_LIBRARY_SOURCES} ${SOURCES} PARENT_SCOPE)
SET (CONTROL_PROBLEM_LIBRARY_HEADERS ${CONTROL_PROBLEM_LIBRARY_HEADERS} ${HEADERS} PARENT_SCOPE)

#  install headers
INSTALL_HEADERS_WITH_DIRECTORY ("${HEADERS}")
