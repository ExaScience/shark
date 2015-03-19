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

# The program to use to edit the cache.
CMAKE_EDIT_COMMAND = /usr/bin/ccmake

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/chakro23/exa2ct/shark/eigen

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/chakro23/exa2ct/shark/eigen/build_dir

# Include any dependencies generated for this target.
include unsupported/test/CMakeFiles/splines.dir/depend.make

# Include the progress variables for this target.
include unsupported/test/CMakeFiles/splines.dir/progress.make

# Include the compile flags for this target's objects.
include unsupported/test/CMakeFiles/splines.dir/flags.make

unsupported/test/CMakeFiles/splines.dir/splines.cpp.o: unsupported/test/CMakeFiles/splines.dir/flags.make
unsupported/test/CMakeFiles/splines.dir/splines.cpp.o: ../unsupported/test/splines.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/chakro23/exa2ct/shark/eigen/build_dir/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object unsupported/test/CMakeFiles/splines.dir/splines.cpp.o"
	cd /home/chakro23/exa2ct/shark/eigen/build_dir/unsupported/test && /opt/gcc-4.8.3/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS)  -DEIGEN_TEST_MAX_SIZE=320 -DEIGEN_TEST_FUNC=splines   -o CMakeFiles/splines.dir/splines.cpp.o -c /home/chakro23/exa2ct/shark/eigen/unsupported/test/splines.cpp

unsupported/test/CMakeFiles/splines.dir/splines.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/splines.dir/splines.cpp.i"
	cd /home/chakro23/exa2ct/shark/eigen/build_dir/unsupported/test && /opt/gcc-4.8.3/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS)  -DEIGEN_TEST_MAX_SIZE=320 -DEIGEN_TEST_FUNC=splines   -E /home/chakro23/exa2ct/shark/eigen/unsupported/test/splines.cpp > CMakeFiles/splines.dir/splines.cpp.i

unsupported/test/CMakeFiles/splines.dir/splines.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/splines.dir/splines.cpp.s"
	cd /home/chakro23/exa2ct/shark/eigen/build_dir/unsupported/test && /opt/gcc-4.8.3/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS)  -DEIGEN_TEST_MAX_SIZE=320 -DEIGEN_TEST_FUNC=splines   -S /home/chakro23/exa2ct/shark/eigen/unsupported/test/splines.cpp -o CMakeFiles/splines.dir/splines.cpp.s

unsupported/test/CMakeFiles/splines.dir/splines.cpp.o.requires:
.PHONY : unsupported/test/CMakeFiles/splines.dir/splines.cpp.o.requires

unsupported/test/CMakeFiles/splines.dir/splines.cpp.o.provides: unsupported/test/CMakeFiles/splines.dir/splines.cpp.o.requires
	$(MAKE) -f unsupported/test/CMakeFiles/splines.dir/build.make unsupported/test/CMakeFiles/splines.dir/splines.cpp.o.provides.build
.PHONY : unsupported/test/CMakeFiles/splines.dir/splines.cpp.o.provides

unsupported/test/CMakeFiles/splines.dir/splines.cpp.o.provides.build: unsupported/test/CMakeFiles/splines.dir/splines.cpp.o

# Object files for target splines
splines_OBJECTS = \
"CMakeFiles/splines.dir/splines.cpp.o"

# External object files for target splines
splines_EXTERNAL_OBJECTS =

unsupported/test/splines: unsupported/test/CMakeFiles/splines.dir/splines.cpp.o
unsupported/test/splines: unsupported/test/CMakeFiles/splines.dir/build.make
unsupported/test/splines: unsupported/test/CMakeFiles/splines.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX executable splines"
	cd /home/chakro23/exa2ct/shark/eigen/build_dir/unsupported/test && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/splines.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
unsupported/test/CMakeFiles/splines.dir/build: unsupported/test/splines
.PHONY : unsupported/test/CMakeFiles/splines.dir/build

unsupported/test/CMakeFiles/splines.dir/requires: unsupported/test/CMakeFiles/splines.dir/splines.cpp.o.requires
.PHONY : unsupported/test/CMakeFiles/splines.dir/requires

unsupported/test/CMakeFiles/splines.dir/clean:
	cd /home/chakro23/exa2ct/shark/eigen/build_dir/unsupported/test && $(CMAKE_COMMAND) -P CMakeFiles/splines.dir/cmake_clean.cmake
.PHONY : unsupported/test/CMakeFiles/splines.dir/clean

unsupported/test/CMakeFiles/splines.dir/depend:
	cd /home/chakro23/exa2ct/shark/eigen/build_dir && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/chakro23/exa2ct/shark/eigen /home/chakro23/exa2ct/shark/eigen/unsupported/test /home/chakro23/exa2ct/shark/eigen/build_dir /home/chakro23/exa2ct/shark/eigen/build_dir/unsupported/test /home/chakro23/exa2ct/shark/eigen/build_dir/unsupported/test/CMakeFiles/splines.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : unsupported/test/CMakeFiles/splines.dir/depend

