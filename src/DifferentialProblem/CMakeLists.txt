# collect header files
SET (HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/AbstractDifferentialProblem.hpp 
             ${CMAKE_CURRENT_SOURCE_DIR}/LinearDifferentialProblem.hpp
             ${CMAKE_CURRENT_SOURCE_DIR}/NonlinearDifferentialProblem.hpp
             ${CMAKE_CURRENT_SOURCE_DIR}/CompositeDifferentialProblem.hpp)

# collect source files
SET (SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/CompositeDifferentialProblem.cpp
             ${CMAKE_CURRENT_SOURCE_DIR}/AbstractDifferentialProblem.cpp)

# create shared library
ADD_LIBRARY (differentialproblem SHARED ${SOURCES} ${HEADERS})

#  install headers
#INSTALL (FILES ${HEADERS} DESTINATION "${HEADER_INSTALL_DIR}")
INSTALL_HEADERS_WITH_DIRECTORY ("${HEADERS}")

# install library
INSTALL (TARGETS differentialproblem DESTINATION ${DYNAMIC_INSTALL_DIR})
