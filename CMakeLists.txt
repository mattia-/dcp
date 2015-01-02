# check on CMake version
CMAKE_MINIMUM_REQUIRED (VERSION 2.6)

# read user defined values for compilers, if present
IF (DEFINED CXX_COMPILER)
    SET (CMAKE_CXX_COMPILER ${CXX_COMPILER})
ENDIF (DEFINED CXX_COMPILER)
IF (DEFINED C_COMPILER)
    SET (CMAKE_C_COMPILER ${C_COMPILER})
ENDIF (DEFINED C_COMPILER)

# set project properties
SET (PROJECT_NAME DCP)
PROJECT (${PROJECT_NAME})

# check on compiler being used
#IF (NOT CMAKE_COMPILER_IS_GNUCXX)
#    MESSAGE (WARNING "-- *** Compiler ID is not GNU. 
#    Proceed, but this may cause problems")
#ENDIF (NOT CMAKE_COMPILER_IS_GNUCXX)
#
#IF (NOT CMAKE_COMPILER_IS_GNUCC)
#    MESSAGE (WARNING "-- *** Linker ID is not GNU. 
#    Proceed, but this may cause problems")
#ENDIF (NOT CMAKE_COMPILER_IS_GNUCC)

# add compilation flags
SET (CMAKE_CXX_FLAGS "-Wall -std=c++0x")

# packages needed
FIND_PACKAGE (DOLFIN REQUIRED)
FIND_PACKAGE (Doxygen)

# add target for Doxygen documentation
IF(DOXYGEN_FOUND)
    CONFIGURE_FILE (${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)
    ADD_CUSTOM_TARGET (doc
        ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        COMMENT "Generating API documentation with Doxygen" VERBATIM
        )
endif(DOXYGEN_FOUND)

# add dolfin specific compiler options
ADD_DEFINITIONS (${DOLFIN_CXX_DEFINITIONS})

# add directories to -I compiler option
INCLUDE_DIRECTORIES (${DOLFIN_INCLUDE_DIRS})
INCLUDE_DIRECTORIES (SYSTEM ${DOLFIN_3RD_PARTY_INCLUDE_DIRS})
INCLUDE_DIRECTORIES (SYSTEM src/)

# no option rpath when building/installing/...
SET (CMAKE_SKIP_RPATH TRUE)

# set directory for installation
IF (NOT DEFINED LIB_INSTALL_DIR)
    SET (LIB_INSTALL_DIR "${CMAKE_SOURCE_DIR}/lib")
ENDIF(NOT DEFINED LIB_INSTALL_DIR)

IF (NOT DEFINED HEADER_INSTALL_DIR)
    SET (HEADER_INSTALL_DIR "${CMAKE_SOURCE_DIR}/include")
ENDIF(NOT DEFINED HEADER_INSTALL_DIR)

# define macro to install header files preserving folder structure
MACRO(INSTALL_HEADERS_WITH_DIRECTORY HEADER_LIST)
    SET (HEADERS ${HEADER_LIST}) # this divides the string passed as input into one variable for each filename
    FOREACH(HEADER ${HEADERS})
        GET_FILENAME_COMPONENT (DIR ${HEADER} PATH) # DIR is now equal to the path of ${HEADER}
                                                    # (that is, we stripped away the name of the file)
        STRING (REGEX REPLACE ".*/src/" "" STRIPPED_DIR ${DIR}) # substitute all the path that leads up to src/ with nothing
        # STRIPPED_DIR now contains the path of the header relative to src/
        INSTALL(FILES ${HEADER} DESTINATION ${HEADER_INSTALL_DIR}/${STRIPPED_DIR})
    ENDFOREACH(HEADER)
ENDMACRO(INSTALL_HEADERS_WITH_DIRECTORY)
# enter subdirectories
ADD_SUBDIRECTORY (src)
