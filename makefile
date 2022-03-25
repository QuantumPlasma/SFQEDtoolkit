#makefile assembled on https://makefiletutorial.com/

MPIC = mpic++
GSL_INST = ${GSL_INSTALLATION_DIR}

CXX = mpic++
#CXXFLAGS = -I/home/smonte/local_installations/openmpi/include
CXXFLAGS = -I${GSL_INST}/include
LDFLAGS = -L${GSL_INST}/lib -lgsl -lgslcblas -lm

OPENMPFLAG = -fopenmp

BUILD_DIR := ./build
TARGET_EXEC := spectrum_tester

SRC_DIRS := ./src
PUSH_DIRS := ./Test_Modules/Pusher_Modules
FUNC_DIRS := ./Test_Modules/Functions

BUILD1 := ./build_push
MAIN1 := ./Test_Modules/Mains/Pusher_Main
TARGET1 := pusher_tester

BUILD2 := ./build_accuracy
MAIN2 := ./Test_Modules/Mains/Accuracy_Main
TARGET2 := accuracy_tester

BUILD3 := ./build_spectrum
MAIN3 := ./Test_Modules/Mains/Spectrum_Main
TARGET3 := spectrum_tester

# Find all the C++ files we want to compile.
# This contains the main SFQEDtoolkit modules.
#SRCS := $(shell find $(SRC_DIRS) -name *.cpp)

SRCS1 := $(shell find $(SRC_DIRS) $(PUSH_DIRS) $(FUNC_DIRS) $(MAIN1) -name *.cpp)
SRCS2 := $(shell find $(SRC_DIRS) $(PUSH_DIRS) $(FUNC_DIRS) $(MAIN2) -name *.cpp)
SRCS3 := $(shell find $(SRC_DIRS) $(PUSH_DIRS) $(FUNC_DIRS) $(MAIN3) -name *.cpp)

# String substitution for every C++ file.
# As an example, hello.cpp turns into ./build/hello.cpp.o
#OBJS := $(SRCS:%=$(BUILD_DIR)/%.o)

OBJS1 := $(SRCS1:%=$(BUILD1)/%.o)
OBJS2 := $(SRCS2:%=$(BUILD2)/%.o)
OBJS3 := $(SRCS3:%=$(BUILD3)/%.o)

# String substitution (suffix version without %).
# As an example, ./build/hello.cpp.o turns into ./build/hello.cpp.d
#(Many build systems add automatically detected make dependencies into the .d file.
# In particular, for C/C++ source files they determine what #include files are
# required and automatically generate that information into the .d file.)
#DEPS := $(OBJS:.o=.d)

DEPS1 := $(OBJS1:.o=.d)
DEPS2 := $(OBJS2:.o=.d)
DEPS3 := $(OBJS3:.o=.d)

# Every folder in ./src will need to be passed to GCC so that it can find header files
INC_DIRS := $(shell find $(SRC_DIRS) $(PUSH_DIRS) $(FUNC_DIRS) -type d)

# Add a prefix to INC_DIRS. So moduleA would become -ImoduleA. GCC understands this -I flag.
# Since we are now adopting a generic folder hierarchy (that is, the headers and source files
# could potentially live in any treelike structure), we need to add all the subdirectories
# within the src folder to the compile path we pass to the compiler (we'll use -I).
INC_FLAGS := $(addprefix -I,$(INC_DIRS))

# The -MMD and -MP flags together generate Makefiles for us!
# These files will have .d instead of .o as the output.
OTHERFLAGS := $(INC_FLAGS) -MMD -MP

#######################################
#first target
pushertester: $(BUILD1)/$(TARGET1)

# The final build step.
$(BUILD1)/$(TARGET1): $(OBJS1)
	@echo Linking and creating executable
	@$(MPIC) $(OPENMPFLAG) $(LDFLAGS) $^ -o $@ #1
	@cp $@ ./$(TARGET1)

# Build step for C++ source
# $(dir $@) is a makefile function that returns
# the target's directory
$(BUILD1)/%.cpp.o: %.cpp
	@mkdir -p $(dir $@)
	@echo Compiling $<
	@$(MPIC) $(CXXFLAGS) $(OPENMPFLAG) $(OTHERFLAGS) -c $< -o $@ #2

#######################################
#second target
accuracytester: $(BUILD2)/$(TARGET2)

# The final build step.
$(BUILD2)/$(TARGET2): $(OBJS2)
	@echo Linking and creating executable
	@$(MPIC) $(OPENMPFLAG) $(LDFLAGS) $^ -o $@ #1
	@cp $@ ./$(TARGET2)

# Build step for C++ source
# $(dir $@) is a makefile function that returns
# the target's directory
$(BUILD2)/%.cpp.o: %.cpp
	@mkdir -p $(dir $@)
	@echo Compiling $<
	@$(MPIC) $(CXXFLAGS) $(OPENMPFLAG) $(OTHERFLAGS) -c $< -o $@ #2

#######################################
#third target
spectrumtester: $(BUILD3)/$(TARGET3)

# The final build step.
$(BUILD3)/$(TARGET3): $(OBJS3)
	@echo Linking and creating executable
	@$(MPIC) $(OPENMPFLAG) $(LDFLAGS) $^ -o $@ #1
	@cp $@ ./$(TARGET3)

# Build step for C++ source
# $(dir $@) is a makefile function that returns
# the target's directory
$(BUILD3)/%.cpp.o: %.cpp
	@mkdir -p $(dir $@)
	@echo Compiling $<
	@$(MPIC) $(CXXFLAGS) $(OPENMPFLAG) $(OTHERFLAGS) -c $< -o $@ #2

.PHONY: clean
clean:
	@echo Cleaning folder
	@rm -r $(BUILD1)
	@rm -r $(BUILD2)
	@rm -r $(BUILD3)

# Include the .d makefiles. The - at the front suppresses the errors of missing
# Makefiles. Initially, all the .d files will be missing, and we don't want those
# errors to show up. Keep in mind that these files are generated by the -MMD -MP
# flags passed to the compiler.
-include $(DEPS)
