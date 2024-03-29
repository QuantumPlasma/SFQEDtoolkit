#makefile assembled on https://makefiletutorial.com/

MPIC = mpic++
CXX = mpic++

FORTRANC = mpif90

#this flag is very important to compile fortran with the c interoperability!!!
LDFLAGS = -lstdc++ #-lgsl -lgslcblas -lm
CPPFLAGS = -std=c++11
NOUNDERSCORE = -fno-underscoring

OPENMPFLAG = -fopenmp

# SRC_DIRS := ./src_SFQEDtoolkit
YOUR_COMPILATION_PATH := ./
FORTRAN_DIRS := ./src_fortran_wrapper

BUILD_FULL_DIR := ./build
BUILD1 := ./build/inner_build
TARGET1 := fortran_tester


# Find all the C++ files we want to compile.
# This contains the main SFQEDtoolkit modules.
SRCS1 := $(shell find $(YOUR_COMPILATION_PATH) -name '*.cpp')
#remember that FORTRAN_DIRS must precede YOUR_COMPILATION_PATH
SRCFORTRAN := $(shell find $(FORTRAN_DIRS) $(YOUR_COMPILATION_PATH) -name '*.f90')


# String substitution for every C++ file.
# As an example, hello.cpp turns into ./build/hello.cpp.o
OBJS1 := $(SRCS1:%=$(BUILD1)/%.o)
OBJSFORTRAN := $(SRCFORTRAN:%=$(BUILD1)/%.o)


# String substitution (suffix version without %).
# As an example, ./build/hello.cpp.o turns into ./build/hello.cpp.d
#(Many build systems add automatically detected make dependencies into the .d file.
# In particular, for C/C++ source files they determine what #include files are
# required and automatically generate that information into the .d file.)
DEPS := $(OBJS:.o=.d)


# Every folder in ./src will need to be passed to GCC so that it can find header files
INC_DIRS := $(shell find $(YOUR_COMPILATION_PATH) -type d)
# INC_DIRS_FORTRAN := $(shell find $(FORTRAN_DIRS) -type d)


# Add a prefix to INC_DIRS. So moduleA would become -ImoduleA. GCC understands this -I flag.
# Since we are now adopting a generic folder hierarchy (that is, the headers and source files
# could potentially live in any treelike structure), we need to add all the subdirectories
# within the src folder to the compile path we pass to the compiler (we'll use -I).
INC_FLAGS := $(addprefix -I,$(INC_DIRS))
# INC_FLAGS_FORTRAN := $(addprefix -I,$(INC_DIRS_FORTRAN))

# The -MMD and -MP flags together generate Makefiles for us!
# These files will have .d instead of .o as the output.
OTHERFLAGS := $(INC_FLAGS) -MMD -MP


#######################################
#first target
fortrantester: $(BUILD1)/$(TARGET1)

# The final build step.
$(BUILD1)/$(TARGET1): $(OBJS1) $(OBJSFORTRAN)
	@echo Linking and creating executable
#	@$(FORTRANC) $(OPENMPFLAG) $(LDFLAGS) $^ -o $@ #1
	@$(FORTRANC) $(LDFLAGS) $^ -o $@ #1
	@cp $@ ./$(TARGET1)

# Build step for C++ source
# $(dir $@) is a makefile function that returns
# the target's directory
$(BUILD1)/%.cpp.o: %.cpp
	@mkdir -p $(dir $@)
	@echo Compiling $<
#	@$(MPIC) $(CXXFLAGS) $(OPENMPFLAG) $(OTHERFLAGS) -c $< -o $@ #2
	@$(MPIC) $(CPPFLAGS) $(OTHERFLAGS) -c $< -o $@ #2

$(BUILD1)/%.f90.o: %.f90
	@mkdir -p $(dir $@)
	@echo Compiling $<
	@$(FORTRANC) $(INC_FLAGS) $(NOUNDERSCORE) -c $< -o $@ #2 

#######################################

.PHONY: clean
clean:
	@echo Cleaning folder
	@rm -r -f $(BUILD_FULL_DIR)
	@rm -f ./*.mod


# Include the .d makefiles. The - at the front suppresses the errors of missing
# Makefiles. Initially, all the .d files will be missing, and we don't want those
# errors to show up. Keep in mind that these files are generated by the -MMD -MP
# flags passed to the compiler.
-include $(DEPS)
