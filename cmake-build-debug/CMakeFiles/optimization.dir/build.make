# CMAKE generated file: DO NOT EDIT!
# Generated by "MinGW Makefiles" Generator, CMake Version 3.23

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

SHELL = cmd.exe

# The CMake executable.
CMAKE_COMMAND = "C:\Program Files\JetBrains\CLion 2022.2.3\bin\cmake\win\bin\cmake.exe"

# The command to remove a file.
RM = "C:\Program Files\JetBrains\CLion 2022.2.3\bin\cmake\win\bin\cmake.exe" -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = C:\Users\szymo\CLionProjects\optimization

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = C:\Users\szymo\CLionProjects\optimization\cmake-build-debug

# Include any dependencies generated for this target.
include CMakeFiles/optimization.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/optimization.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/optimization.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/optimization.dir/flags.make

CMakeFiles/optimization.dir/DRIVER_CODE.cpp.obj: CMakeFiles/optimization.dir/flags.make
CMakeFiles/optimization.dir/DRIVER_CODE.cpp.obj: ../DRIVER_CODE.cpp
CMakeFiles/optimization.dir/DRIVER_CODE.cpp.obj: CMakeFiles/optimization.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=C:\Users\szymo\CLionProjects\optimization\cmake-build-debug\CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/optimization.dir/DRIVER_CODE.cpp.obj"
	C:\PROGRA~1\JETBRA~1\CLION2~1.3\bin\mingw\bin\G__~1.EXE $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/optimization.dir/DRIVER_CODE.cpp.obj -MF CMakeFiles\optimization.dir\DRIVER_CODE.cpp.obj.d -o CMakeFiles\optimization.dir\DRIVER_CODE.cpp.obj -c C:\Users\szymo\CLionProjects\optimization\DRIVER_CODE.cpp

CMakeFiles/optimization.dir/DRIVER_CODE.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/optimization.dir/DRIVER_CODE.cpp.i"
	C:\PROGRA~1\JETBRA~1\CLION2~1.3\bin\mingw\bin\G__~1.EXE $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E C:\Users\szymo\CLionProjects\optimization\DRIVER_CODE.cpp > CMakeFiles\optimization.dir\DRIVER_CODE.cpp.i

CMakeFiles/optimization.dir/DRIVER_CODE.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/optimization.dir/DRIVER_CODE.cpp.s"
	C:\PROGRA~1\JETBRA~1\CLION2~1.3\bin\mingw\bin\G__~1.EXE $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S C:\Users\szymo\CLionProjects\optimization\DRIVER_CODE.cpp -o CMakeFiles\optimization.dir\DRIVER_CODE.cpp.s

CMakeFiles/optimization.dir/classes/matrix.cpp.obj: CMakeFiles/optimization.dir/flags.make
CMakeFiles/optimization.dir/classes/matrix.cpp.obj: ../classes/matrix.cpp
CMakeFiles/optimization.dir/classes/matrix.cpp.obj: CMakeFiles/optimization.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=C:\Users\szymo\CLionProjects\optimization\cmake-build-debug\CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/optimization.dir/classes/matrix.cpp.obj"
	C:\PROGRA~1\JETBRA~1\CLION2~1.3\bin\mingw\bin\G__~1.EXE $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/optimization.dir/classes/matrix.cpp.obj -MF CMakeFiles\optimization.dir\classes\matrix.cpp.obj.d -o CMakeFiles\optimization.dir\classes\matrix.cpp.obj -c C:\Users\szymo\CLionProjects\optimization\classes\matrix.cpp

CMakeFiles/optimization.dir/classes/matrix.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/optimization.dir/classes/matrix.cpp.i"
	C:\PROGRA~1\JETBRA~1\CLION2~1.3\bin\mingw\bin\G__~1.EXE $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E C:\Users\szymo\CLionProjects\optimization\classes\matrix.cpp > CMakeFiles\optimization.dir\classes\matrix.cpp.i

CMakeFiles/optimization.dir/classes/matrix.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/optimization.dir/classes/matrix.cpp.s"
	C:\PROGRA~1\JETBRA~1\CLION2~1.3\bin\mingw\bin\G__~1.EXE $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S C:\Users\szymo\CLionProjects\optimization\classes\matrix.cpp -o CMakeFiles\optimization.dir\classes\matrix.cpp.s

CMakeFiles/optimization.dir/classes/ode_solver.cpp.obj: CMakeFiles/optimization.dir/flags.make
CMakeFiles/optimization.dir/classes/ode_solver.cpp.obj: ../classes/ode_solver.cpp
CMakeFiles/optimization.dir/classes/ode_solver.cpp.obj: CMakeFiles/optimization.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=C:\Users\szymo\CLionProjects\optimization\cmake-build-debug\CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object CMakeFiles/optimization.dir/classes/ode_solver.cpp.obj"
	C:\PROGRA~1\JETBRA~1\CLION2~1.3\bin\mingw\bin\G__~1.EXE $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/optimization.dir/classes/ode_solver.cpp.obj -MF CMakeFiles\optimization.dir\classes\ode_solver.cpp.obj.d -o CMakeFiles\optimization.dir\classes\ode_solver.cpp.obj -c C:\Users\szymo\CLionProjects\optimization\classes\ode_solver.cpp

CMakeFiles/optimization.dir/classes/ode_solver.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/optimization.dir/classes/ode_solver.cpp.i"
	C:\PROGRA~1\JETBRA~1\CLION2~1.3\bin\mingw\bin\G__~1.EXE $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E C:\Users\szymo\CLionProjects\optimization\classes\ode_solver.cpp > CMakeFiles\optimization.dir\classes\ode_solver.cpp.i

CMakeFiles/optimization.dir/classes/ode_solver.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/optimization.dir/classes/ode_solver.cpp.s"
	C:\PROGRA~1\JETBRA~1\CLION2~1.3\bin\mingw\bin\G__~1.EXE $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S C:\Users\szymo\CLionProjects\optimization\classes\ode_solver.cpp -o CMakeFiles\optimization.dir\classes\ode_solver.cpp.s

CMakeFiles/optimization.dir/classes/opt_alg.cpp.obj: CMakeFiles/optimization.dir/flags.make
CMakeFiles/optimization.dir/classes/opt_alg.cpp.obj: ../classes/opt_alg.cpp
CMakeFiles/optimization.dir/classes/opt_alg.cpp.obj: CMakeFiles/optimization.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=C:\Users\szymo\CLionProjects\optimization\cmake-build-debug\CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object CMakeFiles/optimization.dir/classes/opt_alg.cpp.obj"
	C:\PROGRA~1\JETBRA~1\CLION2~1.3\bin\mingw\bin\G__~1.EXE $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/optimization.dir/classes/opt_alg.cpp.obj -MF CMakeFiles\optimization.dir\classes\opt_alg.cpp.obj.d -o CMakeFiles\optimization.dir\classes\opt_alg.cpp.obj -c C:\Users\szymo\CLionProjects\optimization\classes\opt_alg.cpp

CMakeFiles/optimization.dir/classes/opt_alg.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/optimization.dir/classes/opt_alg.cpp.i"
	C:\PROGRA~1\JETBRA~1\CLION2~1.3\bin\mingw\bin\G__~1.EXE $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E C:\Users\szymo\CLionProjects\optimization\classes\opt_alg.cpp > CMakeFiles\optimization.dir\classes\opt_alg.cpp.i

CMakeFiles/optimization.dir/classes/opt_alg.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/optimization.dir/classes/opt_alg.cpp.s"
	C:\PROGRA~1\JETBRA~1\CLION2~1.3\bin\mingw\bin\G__~1.EXE $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S C:\Users\szymo\CLionProjects\optimization\classes\opt_alg.cpp -o CMakeFiles\optimization.dir\classes\opt_alg.cpp.s

CMakeFiles/optimization.dir/classes/solution.cpp.obj: CMakeFiles/optimization.dir/flags.make
CMakeFiles/optimization.dir/classes/solution.cpp.obj: ../classes/solution.cpp
CMakeFiles/optimization.dir/classes/solution.cpp.obj: CMakeFiles/optimization.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=C:\Users\szymo\CLionProjects\optimization\cmake-build-debug\CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building CXX object CMakeFiles/optimization.dir/classes/solution.cpp.obj"
	C:\PROGRA~1\JETBRA~1\CLION2~1.3\bin\mingw\bin\G__~1.EXE $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/optimization.dir/classes/solution.cpp.obj -MF CMakeFiles\optimization.dir\classes\solution.cpp.obj.d -o CMakeFiles\optimization.dir\classes\solution.cpp.obj -c C:\Users\szymo\CLionProjects\optimization\classes\solution.cpp

CMakeFiles/optimization.dir/classes/solution.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/optimization.dir/classes/solution.cpp.i"
	C:\PROGRA~1\JETBRA~1\CLION2~1.3\bin\mingw\bin\G__~1.EXE $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E C:\Users\szymo\CLionProjects\optimization\classes\solution.cpp > CMakeFiles\optimization.dir\classes\solution.cpp.i

CMakeFiles/optimization.dir/classes/solution.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/optimization.dir/classes/solution.cpp.s"
	C:\PROGRA~1\JETBRA~1\CLION2~1.3\bin\mingw\bin\G__~1.EXE $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S C:\Users\szymo\CLionProjects\optimization\classes\solution.cpp -o CMakeFiles\optimization.dir\classes\solution.cpp.s

CMakeFiles/optimization.dir/classes/user_funs.cpp.obj: CMakeFiles/optimization.dir/flags.make
CMakeFiles/optimization.dir/classes/user_funs.cpp.obj: ../classes/user_funs.cpp
CMakeFiles/optimization.dir/classes/user_funs.cpp.obj: CMakeFiles/optimization.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=C:\Users\szymo\CLionProjects\optimization\cmake-build-debug\CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Building CXX object CMakeFiles/optimization.dir/classes/user_funs.cpp.obj"
	C:\PROGRA~1\JETBRA~1\CLION2~1.3\bin\mingw\bin\G__~1.EXE $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/optimization.dir/classes/user_funs.cpp.obj -MF CMakeFiles\optimization.dir\classes\user_funs.cpp.obj.d -o CMakeFiles\optimization.dir\classes\user_funs.cpp.obj -c C:\Users\szymo\CLionProjects\optimization\classes\user_funs.cpp

CMakeFiles/optimization.dir/classes/user_funs.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/optimization.dir/classes/user_funs.cpp.i"
	C:\PROGRA~1\JETBRA~1\CLION2~1.3\bin\mingw\bin\G__~1.EXE $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E C:\Users\szymo\CLionProjects\optimization\classes\user_funs.cpp > CMakeFiles\optimization.dir\classes\user_funs.cpp.i

CMakeFiles/optimization.dir/classes/user_funs.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/optimization.dir/classes/user_funs.cpp.s"
	C:\PROGRA~1\JETBRA~1\CLION2~1.3\bin\mingw\bin\G__~1.EXE $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S C:\Users\szymo\CLionProjects\optimization\classes\user_funs.cpp -o CMakeFiles\optimization.dir\classes\user_funs.cpp.s

# Object files for target optimization
optimization_OBJECTS = \
"CMakeFiles/optimization.dir/DRIVER_CODE.cpp.obj" \
"CMakeFiles/optimization.dir/classes/matrix.cpp.obj" \
"CMakeFiles/optimization.dir/classes/ode_solver.cpp.obj" \
"CMakeFiles/optimization.dir/classes/opt_alg.cpp.obj" \
"CMakeFiles/optimization.dir/classes/solution.cpp.obj" \
"CMakeFiles/optimization.dir/classes/user_funs.cpp.obj"

# External object files for target optimization
optimization_EXTERNAL_OBJECTS =

optimization.exe: CMakeFiles/optimization.dir/DRIVER_CODE.cpp.obj
optimization.exe: CMakeFiles/optimization.dir/classes/matrix.cpp.obj
optimization.exe: CMakeFiles/optimization.dir/classes/ode_solver.cpp.obj
optimization.exe: CMakeFiles/optimization.dir/classes/opt_alg.cpp.obj
optimization.exe: CMakeFiles/optimization.dir/classes/solution.cpp.obj
optimization.exe: CMakeFiles/optimization.dir/classes/user_funs.cpp.obj
optimization.exe: CMakeFiles/optimization.dir/build.make
optimization.exe: CMakeFiles/optimization.dir/linklibs.rsp
optimization.exe: CMakeFiles/optimization.dir/objects1.rsp
optimization.exe: CMakeFiles/optimization.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=C:\Users\szymo\CLionProjects\optimization\cmake-build-debug\CMakeFiles --progress-num=$(CMAKE_PROGRESS_7) "Linking CXX executable optimization.exe"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles\optimization.dir\link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/optimization.dir/build: optimization.exe
.PHONY : CMakeFiles/optimization.dir/build

CMakeFiles/optimization.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles\optimization.dir\cmake_clean.cmake
.PHONY : CMakeFiles/optimization.dir/clean

CMakeFiles/optimization.dir/depend:
	$(CMAKE_COMMAND) -E cmake_depends "MinGW Makefiles" C:\Users\szymo\CLionProjects\optimization C:\Users\szymo\CLionProjects\optimization C:\Users\szymo\CLionProjects\optimization\cmake-build-debug C:\Users\szymo\CLionProjects\optimization\cmake-build-debug C:\Users\szymo\CLionProjects\optimization\cmake-build-debug\CMakeFiles\optimization.dir\DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/optimization.dir/depend
