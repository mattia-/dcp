# collect header files
SET (HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/AbstractOptimizer.hpp
             ${CMAKE_CURRENT_SOURCE_DIR}/BacktrackingOptimizer.hpp
             ${CMAKE_CURRENT_SOURCE_DIR}/DirichletControlValueUpdater.hpp
             ${CMAKE_CURRENT_SOURCE_DIR}/DistributedControlValueUpdater.hpp
             ${CMAKE_CURRENT_SOURCE_DIR}/NeumannControlValueUpdater.hpp
             ${CMAKE_CURRENT_SOURCE_DIR}/Optimizer.hpp
             ${CMAKE_CURRENT_SOURCE_DIR}/dotproduct.h)

# collect source files
SET (SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/AbstractOptimizer.cpp
             ${CMAKE_CURRENT_SOURCE_DIR}/BacktrackingOptimizer.cpp
             ${CMAKE_CURRENT_SOURCE_DIR}/DirichletControlValueUpdater.cpp
             ${CMAKE_CURRENT_SOURCE_DIR}/DistributedControlValueUpdater.cpp
             ${CMAKE_CURRENT_SOURCE_DIR}/NeumannControlValueUpdater.cpp)

# add files to shared library dependancies
SET (CONTROL_PROBLEM_LIBRARY_SOURCES ${CONTROL_PROBLEM_LIBRARY_SOURCES} ${SOURCES} PARENT_SCOPE)
SET (CONTROL_PROBLEM_LIBRARY_HEADERS ${CONTROL_PROBLEM_LIBRARY_HEADERS} ${HEADERS} PARENT_SCOPE)

#  install headers
INSTALL_HEADERS_WITH_DIRECTORY ("${HEADERS}")
