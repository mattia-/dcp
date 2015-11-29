# collect header files
SET (HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/GenericProblem.h 
             ${CMAKE_CURRENT_SOURCE_DIR}/GenericLinearProblem.h
             ${CMAKE_CURRENT_SOURCE_DIR}/GenericNonlinearProblem.h
             ${CMAKE_CURRENT_SOURCE_DIR}/LinearProblem.h
             ${CMAKE_CURRENT_SOURCE_DIR}/NonlinearProblem.h
             ${CMAKE_CURRENT_SOURCE_DIR}/TimeDependentProblem.h
             ${CMAKE_CURRENT_SOURCE_DIR}/AlgebraicProblem.h
             ${CMAKE_CURRENT_SOURCE_DIR}/GenericEquationSystem.h
             ${CMAKE_CURRENT_SOURCE_DIR}/EquationSystem.h
             ${CMAKE_CURRENT_SOURCE_DIR}/TimeDependentEquationSystem.h
             ${CMAKE_CURRENT_SOURCE_DIR}/SubdomainType.h
             ${CMAKE_CURRENT_SOURCE_DIR}/problems.h)

# collect source files
SET (SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/EquationSystem.cpp
             ${CMAKE_CURRENT_SOURCE_DIR}/TimeDependentEquationSystem.cpp
             ${CMAKE_CURRENT_SOURCE_DIR}/TimeDependentProblem.cpp
             ${CMAKE_CURRENT_SOURCE_DIR}/GenericProblem.cpp
             ${CMAKE_CURRENT_SOURCE_DIR}/GenericLinearProblem.cpp
             ${CMAKE_CURRENT_SOURCE_DIR}/GenericNonlinearProblem.cpp
             ${CMAKE_CURRENT_SOURCE_DIR}/AlgebraicProblem.cpp
             ${CMAKE_CURRENT_SOURCE_DIR}/GenericEquationSystem.cpp)

# add files to shared library dependancies
SET (CONTROL_PROBLEM_LIBRARY_SOURCES ${CONTROL_PROBLEM_LIBRARY_SOURCES} ${SOURCES} PARENT_SCOPE)
SET (CONTROL_PROBLEM_LIBRARY_HEADERS ${CONTROL_PROBLEM_LIBRARY_HEADERS} ${HEADERS} PARENT_SCOPE)

#  install headers
INSTALL_HEADERS_WITH_DIRECTORY ("${HEADERS}")
