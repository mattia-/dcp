# collect header files
SET (HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/VariableExpression.hpp)

# collect source files
SET (SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/VariableExpression.cpp)

# create shared library
ADD_LIBRARY (functional SHARED ${SOURCES} ${HEADERS})

#  install headers
#INSTALL (FILES ${HEADERS} DESTINATION "${HEADER_INSTALL_DIR}")
INSTALL_HEADERS_WITH_DIRECTORY ("${HEADERS}")

# install library
INSTALL (TARGETS functional DESTINATION ${DYNAMIC_INSTALL_DIR})
