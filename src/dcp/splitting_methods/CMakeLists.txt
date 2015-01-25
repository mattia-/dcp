# collect header files
SET (HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/AbstractSplittingMethod.h 
             ${CMAKE_CURRENT_SOURCE_DIR}/NavierStokesSplittingMethod.h
             ${CMAKE_CURRENT_SOURCE_DIR}/ChorinTemamMethod.h
             ${CMAKE_CURRENT_SOURCE_DIR}/IncrementalChorinTemamMethod.h
             ${CMAKE_CURRENT_SOURCE_DIR}/splitting_methods.h)

# collect source files
SET (SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/AbstractSplittingMethod.cpp
             ${CMAKE_CURRENT_SOURCE_DIR}/NavierStokesSplittingMethod.cpp)

# add files to shared library dependancies
SET (CONTROL_PROBLEM_LIBRARY_SOURCES ${CONTROL_PROBLEM_LIBRARY_SOURCES} ${SOURCES} PARENT_SCOPE)
SET (CONTROL_PROBLEM_LIBRARY_HEADERS ${CONTROL_PROBLEM_LIBRARY_HEADERS} ${HEADERS} PARENT_SCOPE)

#  install headers
INSTALL_HEADERS_WITH_DIRECTORY ("${HEADERS}")
