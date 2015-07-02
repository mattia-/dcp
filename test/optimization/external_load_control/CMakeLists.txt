SET (PROJECT_NAME optimization_external_load_control)
PROJECT (${PROJECT_NAME})

SET (CMAKE_SKIP_RPATH TRUE)

FILE (READ "../../cmake_header_dir.install" INCLUDE_DIRS)
FILE (READ "../../cmake_lib_dir.install" LINK_DIRS)

ADD_SUBDIRECTORY (src)
