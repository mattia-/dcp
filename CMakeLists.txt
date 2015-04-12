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
SET (PROJECT_NAME dcp)
PROJECT (${PROJECT_NAME})

# get packages needed by dcp
# dolfin
FIND_PACKAGE (DOLFIN REQUIRED)
INCLUDE (${DOLFIN_USE_FILE})
# doxygen
FIND_PACKAGE (Doxygen)

# check if file "./src/dcp/optimizers" exists and was compiler for
# the correct version of dolfin. Remember that the version of dolfin
# is defined by the dolfin cmake package and is stored in the variable
# DOLFIN_VERSION_STR
# </====================================
SET (DOT_PRODUCT_FILE "${CMAKE_CURRENT_SOURCE_DIR}/src/dcp/optimizers/dotproduct.h")
# if the file does not exists, compile it
IF (NOT EXISTS ${DOT_PRODUCT_FILE})
    MESSAGE (STATUS "Compiling ufl file \"dotproduct.ufl\"")
    EXECUTE_PROCESS (COMMAND ffc -l dolfin -O dotproduct.ufl
                     WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/src/dcp/optimizers/
                     OUTPUT_QUIET)
# if the file exists, check the versions
ELSE (NOT EXISTS ${DOT_PRODUCT_FILE})
    MESSAGE (STATUS "Checking compatibility of file \"${DOT_PRODUCT_FILE}\"")

    # get FFC version from file
    FILE (STRINGS ${DOT_PRODUCT_FILE} DOT_PRODUCT_FILE_FFC_VERSION
          REGEX "FFC version [0-9.]*$")
    STRING (REGEX MATCH "[0-9.]*$" DOT_PRODUCT_FFC_VERSION ${DOT_PRODUCT_FILE_FFC_VERSION})

    # compare dolfin version with ffc version extracted from file.
    # We need to add a dot at the end of ${DOLFIN_VERSION_STR} since the regexp
    # in the last expression matched also the dot at the end of the ffc version in
    # "dotproduct.h"
    STRING (COMPARE EQUAL "${DOLFIN_VERSION_STR}." ${DOT_PRODUCT_FFC_VERSION} VERSION_COMPARISON)

    # if the output is 0, then the versions do not match. In this case, recompile
    IF (${VERSION_COMPARISON} EQUAL 0)
        MESSAGE (STATUS "Mismatch between system FFC version and \"dotproduct.h\" FFC version")
        MESSAGE (STATUS "Compiling ufl file \"dotproduct.ufl\"")
        EXECUTE_PROCESS (COMMAND ffc -l dolfin -O dotproduct.ufl
                         WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/src/dcp/optimizers/
                         OUTPUT_QUIET)
    ENDIF (${VERSION_COMPARISON} EQUAL 0)

    MESSAGE (STATUS "Checking compatibility of file \"${DOT_PRODUCT_FILE}\" -- done")
ENDIF (NOT EXISTS ${DOT_PRODUCT_FILE})
# ====================================/>



# add compilation flags
SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")

# add target for Doxygen documentation
# </====================================
IF (DOXYGEN_FOUND)
    CONFIGURE_FILE (${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)
    ADD_CUSTOM_TARGET (doc
        ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        COMMENT "Generating API documentation with Doxygen" VERBATIM)
ENDIF (DOXYGEN_FOUND)
# ====================================/>

# add target for uninstallation
# </====================================
CONFIGURE_FILE ("${CMAKE_CURRENT_SOURCE_DIR}/cmake_uninstall.cmake.in"
                "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake"
                IMMEDIATE @ONLY)
ADD_CUSTOM_TARGET (uninstall
                   COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake)
# ====================================/>

# add dcp specific directories to -I compiler option
INCLUDE_DIRECTORIES (SYSTEM src/)

# no option rpath when building/installing/...
SET (CMAKE_SKIP_RPATH TRUE)

# set directory for installation
IF (NOT DEFINED LIB_INSTALL_DIR)
    SET (LIB_INSTALL_DIR "${CMAKE_INSTALL_PREFIX}/lib")
ENDIF (NOT DEFINED LIB_INSTALL_DIR)

IF (NOT DEFINED HEADER_INSTALL_DIR)
    SET (HEADER_INSTALL_DIR "${CMAKE_INSTALL_PREFIX}/include")
ENDIF (NOT DEFINED HEADER_INSTALL_DIR)

MESSAGE (STATUS "Library installation directory is ${LIB_INSTALL_DIR}")
MESSAGE (STATUS "Header files installation directory is ${HEADER_INSTALL_DIR}")

# define macro to install header files preserving folder structure
# </====================================
MACRO(INSTALL_HEADERS_WITH_DIRECTORY HEADER_LIST)
    SET (HEADERS ${HEADER_LIST}) # this divides the string passed as input into one variable for each filename
    FOREACH(HEADER ${HEADERS})
        GET_FILENAME_COMPONENT (DIR ${HEADER} PATH) # DIR is now equal to the path of ${HEADER}
                                                    # (that is, we stripped away the name of the file)
        STRING (REGEX REPLACE ".*/src/" "" STRIPPED_DIR ${DIR}) # substitute all the path that leads up to src/ with nothing
        # STRIPPED_DIR now contains the path of the header relative to src/
        INSTALL(FILES ${HEADER} DESTINATION ${HEADER_INSTALL_DIR}/${STRIPPED_DIR})
    ENDFOREACH(HEADER)
ENDMACRO (INSTALL_HEADERS_WITH_DIRECTORY)
# ====================================/>

# enter subdirectories
ADD_SUBDIRECTORY (src)
