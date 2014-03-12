# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 2.8

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list

# Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/mattia/Dropbox/polimi/tesi/progetto_PACS/git_repo/src

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/mattia/Dropbox/polimi/tesi/progetto_PACS/git_repo/src/build

# Include any dependencies generated for this target.
include Factory/CMakeFiles/linearsolverproxy.dir/depend.make

# Include the progress variables for this target.
include Factory/CMakeFiles/linearsolverproxy.dir/progress.make

# Include the compile flags for this target's objects.
include Factory/CMakeFiles/linearsolverproxy.dir/flags.make

Factory/CMakeFiles/linearsolverproxy.dir/LinearSolverProxy.cpp.o: Factory/CMakeFiles/linearsolverproxy.dir/flags.make
Factory/CMakeFiles/linearsolverproxy.dir/LinearSolverProxy.cpp.o: ../Factory/LinearSolverProxy.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/mattia/Dropbox/polimi/tesi/progetto_PACS/git_repo/src/build/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object Factory/CMakeFiles/linearsolverproxy.dir/LinearSolverProxy.cpp.o"
	cd /home/mattia/Dropbox/polimi/tesi/progetto_PACS/git_repo/src/build/Factory && /usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/linearsolverproxy.dir/LinearSolverProxy.cpp.o -c /home/mattia/Dropbox/polimi/tesi/progetto_PACS/git_repo/src/Factory/LinearSolverProxy.cpp

Factory/CMakeFiles/linearsolverproxy.dir/LinearSolverProxy.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/linearsolverproxy.dir/LinearSolverProxy.cpp.i"
	cd /home/mattia/Dropbox/polimi/tesi/progetto_PACS/git_repo/src/build/Factory && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/mattia/Dropbox/polimi/tesi/progetto_PACS/git_repo/src/Factory/LinearSolverProxy.cpp > CMakeFiles/linearsolverproxy.dir/LinearSolverProxy.cpp.i

Factory/CMakeFiles/linearsolverproxy.dir/LinearSolverProxy.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/linearsolverproxy.dir/LinearSolverProxy.cpp.s"
	cd /home/mattia/Dropbox/polimi/tesi/progetto_PACS/git_repo/src/build/Factory && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/mattia/Dropbox/polimi/tesi/progetto_PACS/git_repo/src/Factory/LinearSolverProxy.cpp -o CMakeFiles/linearsolverproxy.dir/LinearSolverProxy.cpp.s

Factory/CMakeFiles/linearsolverproxy.dir/LinearSolverProxy.cpp.o.requires:
.PHONY : Factory/CMakeFiles/linearsolverproxy.dir/LinearSolverProxy.cpp.o.requires

Factory/CMakeFiles/linearsolverproxy.dir/LinearSolverProxy.cpp.o.provides: Factory/CMakeFiles/linearsolverproxy.dir/LinearSolverProxy.cpp.o.requires
	$(MAKE) -f Factory/CMakeFiles/linearsolverproxy.dir/build.make Factory/CMakeFiles/linearsolverproxy.dir/LinearSolverProxy.cpp.o.provides.build
.PHONY : Factory/CMakeFiles/linearsolverproxy.dir/LinearSolverProxy.cpp.o.provides

Factory/CMakeFiles/linearsolverproxy.dir/LinearSolverProxy.cpp.o.provides.build: Factory/CMakeFiles/linearsolverproxy.dir/LinearSolverProxy.cpp.o

# Object files for target linearsolverproxy
linearsolverproxy_OBJECTS = \
"CMakeFiles/linearsolverproxy.dir/LinearSolverProxy.cpp.o"

# External object files for target linearsolverproxy
linearsolverproxy_EXTERNAL_OBJECTS =

Factory/liblinearsolverproxy.so: Factory/CMakeFiles/linearsolverproxy.dir/LinearSolverProxy.cpp.o
Factory/liblinearsolverproxy.so: Factory/CMakeFiles/linearsolverproxy.dir/build.make
Factory/liblinearsolverproxy.so: Factory/CMakeFiles/linearsolverproxy.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX shared library liblinearsolverproxy.so"
	cd /home/mattia/Dropbox/polimi/tesi/progetto_PACS/git_repo/src/build/Factory && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/linearsolverproxy.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
Factory/CMakeFiles/linearsolverproxy.dir/build: Factory/liblinearsolverproxy.so
.PHONY : Factory/CMakeFiles/linearsolverproxy.dir/build

Factory/CMakeFiles/linearsolverproxy.dir/requires: Factory/CMakeFiles/linearsolverproxy.dir/LinearSolverProxy.cpp.o.requires
.PHONY : Factory/CMakeFiles/linearsolverproxy.dir/requires

Factory/CMakeFiles/linearsolverproxy.dir/clean:
	cd /home/mattia/Dropbox/polimi/tesi/progetto_PACS/git_repo/src/build/Factory && $(CMAKE_COMMAND) -P CMakeFiles/linearsolverproxy.dir/cmake_clean.cmake
.PHONY : Factory/CMakeFiles/linearsolverproxy.dir/clean

Factory/CMakeFiles/linearsolverproxy.dir/depend:
	cd /home/mattia/Dropbox/polimi/tesi/progetto_PACS/git_repo/src/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/mattia/Dropbox/polimi/tesi/progetto_PACS/git_repo/src /home/mattia/Dropbox/polimi/tesi/progetto_PACS/git_repo/src/Factory /home/mattia/Dropbox/polimi/tesi/progetto_PACS/git_repo/src/build /home/mattia/Dropbox/polimi/tesi/progetto_PACS/git_repo/src/build/Factory /home/mattia/Dropbox/polimi/tesi/progetto_PACS/git_repo/src/build/Factory/CMakeFiles/linearsolverproxy.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : Factory/CMakeFiles/linearsolverproxy.dir/depend

