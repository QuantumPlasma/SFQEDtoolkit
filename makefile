MPIC = mpic++
GSL_INST = /home/smonte/extend/gsl_installation_dir
MPI_INST = /home/smonte/local_installations/openmpi

CXX = mpic++
#CXXFLAGS = -I/home/smonte/local_installations/openmpi/include
CXXFLAGS = -I/home/smonte/extend/gsl_installation_dir/include
LDFLAGS = -L/home/smonte/extend/gsl_installation_dir/lib -lgsl -lgslcblas -lm

OPENMPFLAG = -fopenmp
MAINFILEOBJ = main.o
CPPFILEOBJS = Chebyshev_1D.o Test_Functions.o Chebyshev_Utils.o

BUILD_DIR := ./build
SRC_DIRS := ./src
TARGET_EXEC := spectrum_tester

#SFQEDUSER_DIRS = /home/smonte/extend/sfqed_tk_user/include
#SFQEDUSERFLAGS := $(addprefix -I,$(SFQEDUSER_DIRS))

# Find all the C++ files we want to compile
SRCS := $(shell find $(SRC_DIRS) -name *.cpp)

# String substitution for every C++ file.
# As an example, hello.cpp turns into ./build/hello.cpp.o
OBJS := $(SRCS:%=$(BUILD_DIR)/%.o)

# String substitution (suffix version without %).
# As an example, ./build/hello.cpp.o turns into ./build/hello.cpp.d
#(Many build systems add automatically detected make dependencies into the .d file.
# In particular, for C/C++ source files they determine what #include files are
# required and automatically generate that information into the .d file.)
DEPS := $(OBJS:.o=.d)

# Every folder in ./src will need to be passed to GCC so that it can find header files
INC_DIRS := $(shell find $(SRC_DIRS) -type d)
# Add a prefix to INC_DIRS. So moduleA would become -ImoduleA. GCC understands this -I flag.
# Since we are now adopting a generic folder hierarchy (that is, the headers and source files
# could potentially live in any treelike structure), we need to add all the subdirectories
# within the src folder to the compile path we pass to the compiler (we'll use -I).
INC_FLAGS := $(addprefix -I,$(INC_DIRS))

# The -MMD and -MP flags together generate Makefiles for us!
# These files will have .d instead of .o as the output.
OTHERFLAGS := $(INC_FLAGS) -MMD -MP

# The final build step.
$(BUILD_DIR)/$(TARGET_EXEC): $(OBJS)
	@echo Linking and creating executable
	@$(MPIC) $(OPENMPFLAG) $(LDFLAGS) $^ -o $@ #1
	@cp $@ ./$(TARGET_EXEC)

# Build step for C++ source
# $(dir $@) is a makefile function that returns
# the target's directory
$(BUILD_DIR)/%.cpp.o: %.cpp
	@mkdir -p $(dir $@)
	@echo Compiling $<
	@$(MPIC) $(CXXFLAGS) $(OPENMPFLAG) $(OTHERFLAGS) -c $< -o $@ #2

.PHONY: clean
clean:
	@echo Cleaning folder
	@rm -r $(BUILD_DIR)

# Include the .d makefiles. The - at the front suppresses the errors of missing
# Makefiles. Initially, all the .d files will be missing, and we don't want those
# errors to show up.
-include $(DEPS)
