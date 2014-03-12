# Install script for directory: /home/mattia/Dropbox/polimi/tesi/progetto_PACS/git_repo/src/Factory

# Set the install prefix
IF(NOT DEFINED CMAKE_INSTALL_PREFIX)
  SET(CMAKE_INSTALL_PREFIX "/usr/local")
ENDIF(NOT DEFINED CMAKE_INSTALL_PREFIX)
STRING(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
IF(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  IF(BUILD_TYPE)
    STRING(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  ELSE(BUILD_TYPE)
    SET(CMAKE_INSTALL_CONFIG_NAME "")
  ENDIF(BUILD_TYPE)
  MESSAGE(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
ENDIF(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)

# Set the component getting installed.
IF(NOT CMAKE_INSTALL_COMPONENT)
  IF(COMPONENT)
    MESSAGE(STATUS "Install component: \"${COMPONENT}\"")
    SET(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  ELSE(COMPONENT)
    SET(CMAKE_INSTALL_COMPONENT)
  ENDIF(COMPONENT)
ENDIF(NOT CMAKE_INSTALL_COMPONENT)

# Install shared libraries without execute permission?
IF(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  SET(CMAKE_INSTALL_SO_NO_EXE "1")
ENDIF(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  list(APPEND CPACK_ABSOLUTE_DESTINATION_FILES
   "/home/mattia/Dropbox/polimi/tesi/progetto_PACS/git_repo/src/../lib/liblinearsolverproxy.so")
FILE(INSTALL DESTINATION "/home/mattia/Dropbox/polimi/tesi/progetto_PACS/git_repo/src/../lib" TYPE SHARED_LIBRARY FILES "/home/mattia/Dropbox/polimi/tesi/progetto_PACS/git_repo/src/build/Factory/liblinearsolverproxy.so")
  IF(EXISTS "$ENV{DESTDIR}/home/mattia/Dropbox/polimi/tesi/progetto_PACS/git_repo/src/../lib/liblinearsolverproxy.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/home/mattia/Dropbox/polimi/tesi/progetto_PACS/git_repo/src/../lib/liblinearsolverproxy.so")
    IF(CMAKE_INSTALL_DO_STRIP)
      EXECUTE_PROCESS(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}/home/mattia/Dropbox/polimi/tesi/progetto_PACS/git_repo/src/../lib/liblinearsolverproxy.so")
    ENDIF(CMAKE_INSTALL_DO_STRIP)
  ENDIF()
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  list(APPEND CPACK_ABSOLUTE_DESTINATION_FILES
   "/home/mattia/Dropbox/polimi/tesi/progetto_PACS/git_repo/src/../include/GenericFactory.hpp;/home/mattia/Dropbox/polimi/tesi/progetto_PACS/git_repo/src/../include/LinearSolverFactory.hpp;/home/mattia/Dropbox/polimi/tesi/progetto_PACS/git_repo/src/../include/Proxy.hpp")
FILE(INSTALL DESTINATION "/home/mattia/Dropbox/polimi/tesi/progetto_PACS/git_repo/src/../include" TYPE FILE FILES
    "/home/mattia/Dropbox/polimi/tesi/progetto_PACS/git_repo/src/Factory/GenericFactory.hpp"
    "/home/mattia/Dropbox/polimi/tesi/progetto_PACS/git_repo/src/Factory/LinearSolverFactory.hpp"
    "/home/mattia/Dropbox/polimi/tesi/progetto_PACS/git_repo/src/Factory/Proxy.hpp"
    )
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

