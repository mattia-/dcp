# collect header files
SET (HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/AbstractProblem.h 
             ${CMAKE_CURRENT_SOURCE_DIR}/LinearProblem.h
             ${CMAKE_CURRENT_SOURCE_DIR}/MeshManager.h
             ${CMAKE_CURRENT_SOURCE_DIR}/MovingAbstractProblem.h
             ${CMAKE_CURRENT_SOURCE_DIR}/AdditionalProcessor.h
             ${CMAKE_CURRENT_SOURCE_DIR}/MovingLinearProblem.h
             ${CMAKE_CURRENT_SOURCE_DIR}/MovingTimeDependentProblem.h
             ${CMAKE_CURRENT_SOURCE_DIR}/DefaultPostProcessor.h
             ${CMAKE_CURRENT_SOURCE_DIR}/NonlinearProblem.h
             ${CMAKE_CURRENT_SOURCE_DIR}/TimeDependentProblem.h
             ${CMAKE_CURRENT_SOURCE_DIR}/AbstractEquationSystem.h
             ${CMAKE_CURRENT_SOURCE_DIR}/EquationSystem.h
             ${CMAKE_CURRENT_SOURCE_DIR}/TimeDependentEquationSystem.h
             ${CMAKE_CURRENT_SOURCE_DIR}/SubdomainType.h
             ${CMAKE_CURRENT_SOURCE_DIR}/utilities.h
             ${CMAKE_CURRENT_SOURCE_DIR}/meshManagerUflTools.h
             ${CMAKE_CURRENT_SOURCE_DIR}/GetPot.h
             ${CMAKE_CURRENT_SOURCE_DIR}/differential_problems.h)

# collect source files
SET (SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/EquationSystem.cpp
             ${CMAKE_CURRENT_SOURCE_DIR}/TimeDependentEquationSystem.cpp
             ${CMAKE_CURRENT_SOURCE_DIR}/TimeDependentProblem.cpp
             ${CMAKE_CURRENT_SOURCE_DIR}/AbstractProblem.cpp
             ${CMAKE_CURRENT_SOURCE_DIR}/MovingAbstractProblem.cpp
             ${CMAKE_CURRENT_SOURCE_DIR}/MovingTimeDependentProblem.cpp
             ${CMAKE_CURRENT_SOURCE_DIR}/utilities.cpp
             ${CMAKE_CURRENT_SOURCE_DIR}/AbstractEquationSystem.cpp)

# add files to shared library dependancies
SET (CONTROL_PROBLEM_LIBRARY_SOURCES ${CONTROL_PROBLEM_LIBRARY_SOURCES} ${SOURCES} PARENT_SCOPE)
SET (CONTROL_PROBLEM_LIBRARY_HEADERS ${CONTROL_PROBLEM_LIBRARY_HEADERS} ${HEADERS} PARENT_SCOPE)

#  install headers
INSTALL_HEADERS_WITH_DIRECTORY ("${HEADERS}")
