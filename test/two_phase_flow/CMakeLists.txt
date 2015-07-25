SET (PROJECT_NAME time_dependent_cavity_NS)
PROJECT (${PROJECT_NAME})

SET (CMAKE_SKIP_RPATH TRUE)

FILE (READ "../../../cmake_header_dir.install" INCLUDE_DIRS)
FILE (READ "../../../cmake_lib_dir.install" LINK_DIRS)

ADD_SUBDIRECTORY (src)
