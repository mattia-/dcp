# collect header files
SET (HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/SubdomainType.hpp
             ${CMAKE_CURRENT_SOURCE_DIR}/ComparisonFunctors.hpp) 

# collect source files
SET (SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/ComparisonFunctors.cpp)

# create static library
ADD_LIBRARY (comparisonfunctors ${SOURCES} ${HEADERS})

#  install headers
#INSTALL (FILES ${HEADERS} DESTINATION "${HEADER_INSTALL_DIR}")
INSTALL_HEADERS_WITH_DIRECTORY ("${HEADERS}")

# install library
INSTALL (TARGETS comparisonfunctors DESTINATION ${STATIC_INSTALL_DIR})
