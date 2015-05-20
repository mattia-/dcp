# collect header files
SET (HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/AbstractDescentMethod.h
             ${CMAKE_CURRENT_SOURCE_DIR}/BacktrackingOptimizer.h
             ${CMAKE_CURRENT_SOURCE_DIR}/DirichletControlUpdater.h
             ${CMAKE_CURRENT_SOURCE_DIR}/DistributedControlUpdater.h
             ${CMAKE_CURRENT_SOURCE_DIR}/NeumannControlUpdater.h
             ${CMAKE_CURRENT_SOURCE_DIR}/GradientSearchDirection.h
             ${CMAKE_CURRENT_SOURCE_DIR}/optimizers.h)

# collect source files
SET (SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/AbstractDescentMethod.cpp
             ${CMAKE_CURRENT_SOURCE_DIR}/BacktrackingOptimizer.cpp
             ${CMAKE_CURRENT_SOURCE_DIR}/DirichletControlUpdater.cpp
             ${CMAKE_CURRENT_SOURCE_DIR}/DistributedControlUpdater.cpp
             ${CMAKE_CURRENT_SOURCE_DIR}/NeumannControlUpdater.cpp
             ${CMAKE_CURRENT_SOURCE_DIR}/GradientSearchDirection.cpp)

# add files to shared library dependancies
SET (CONTROL_PROBLEM_LIBRARY_SOURCES ${CONTROL_PROBLEM_LIBRARY_SOURCES} ${SOURCES} PARENT_SCOPE)
SET (CONTROL_PROBLEM_LIBRARY_HEADERS ${CONTROL_PROBLEM_LIBRARY_HEADERS} ${HEADERS} PARENT_SCOPE)

#  install headers
INSTALL_HEADERS_WITH_DIRECTORY ("${HEADERS}")
