#=========================================================================
# Gamma top level Makefile
#=========================================================================

include Makefile.config

SRCS = 	arr.cpp\
	AudioIO.cpp\
	Conversion.cpp\
	DFT.cpp\
	FFT_fftpack.cpp\
	fftpack++1.cpp\
	fftpack++2.cpp\
	File.cpp\
	Print.cpp\
	scl.cpp\
	SoundFile.cpp\
	Sync.cpp

OBJS = $(SRCS:.cpp=.o)
OBJS := $(addprefix $(OBJ_DIR), $(OBJS))
SRCS := $(addprefix $(SRC_DIR), $(SRCS))

CPPFLAGS	+= $(EXT_CPPFLAGS)
LDFLAGS		+= $(EXT_LDFLAGS)

CPPFLAGS	+= $(addprefix -I, $(INC_DIRS) $(RINC_DIRS))
LDFLAGS		:= $(addprefix -L, $(MY_LIB_DIRS)) $(LDFLAGS)

LDFLAGS		+= -L$(EXT_LIB_DIR)

DEPS		:= $(OBJS:.o=.d)
CFLAGS		:= $(CPPFLAGS) $(CFLAGS) $(CXXFLAGS)

#--------------------------------------------------------------------------
# Rules
#--------------------------------------------------------------------------

# Force these targets to always execute
.PHONY: clean cleanall external tests

# Build object file from C++ source
$(OBJ_DIR)%.o: %.cpp
	@echo CXX $< $@
	@$(CXX) -c $(CFLAGS) $< -o $@

# Create dependency file from C++ source
$(OBJ_DIR)%.d: %.cpp | createFolders
	@set -e; rm -f $@; \
	$(CXX) -MM $(CPPFLAGS) $< > $@.$$$$; \
	sed 's,\($*\)\.o[ :]*,$(OBJ_DIR)\1.o : ,g' < $@.$$$$ > $@; \
	rm -f $@.$$$$
#	sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $@.$$$$ > $@;

# Build static library
$(SLIB_PATH): createFolders $(OBJS)
	@echo AR $@
	@$(AR) $@ $(OBJS)
	@$(RANLIB) $@

# Build dynamic (shared) library
$(DLIB_FILE): createFolders external $(OBJS)
	@echo SH $@
	@$(CXX) $(DLIBFLAGS) $(LDFLAGS) -o $@ $(OBJS)

# Create file with settings for linking to external libraries
external:
	@echo '\
CPPFLAGS += $(EXT_CPPFLAGS) \r\n\
LDFLAGS  += $(EXT_LDFLAGS) \
'> Makefile.external

# Compile and run source files in examples/ and tests/ folders
examples/%.cpp tests/%.cpp: $(SLIB_PATH) FORCE
	@$(CXX) $(CFLAGS) -o $(BIN_DIR)$(*F) $@ $(LDFLAGS) $(SLIB_PATH)
ifneq ($(AUTORUN), 0)
	@$(BIN_DIR)$(*F)
endif

tests: $(SLIB_PATH)
	@$(MAKE) -C $(TEST_DIR)

# Remove active build configuration binary files
clean:
	@$(RM) $(OBJ_DIR)* $(OBJ_DIR) $(BIN_DIR)* $(BIN_DIR)

# Remove all build configuration binary files
cleanall:
	@$(MAKE) clean BUILD_CONFIG=release
	@$(MAKE) clean BUILD_CONFIG=debug
	@$(RM) $(BUILD_DIR)* $(BUILD_DIR)

# Clean and rebuild library
rebuild: clean $(SLIB_PATH)

# Install library into path specified by INSTALL_DIR
# Include files are copied into INSTALL_DIR/include/LIB_NAME and
# library files are copied to INSTALL_DIR/lib
install: $(SLIB_PATH)
#	@echo 'INSTALL $(INSTALL_DIR)'
	@$(INSTALL) -d $(INSTALL_DIR)
	@$(INSTALL) -d $(INSTALL_DIR)lib
	@$(INSTALL) -d $(INSTALL_DIR)include/$(LIB_NAME)
	@$(INSTALL) -c -m 644 $(SLIB_PATH) $(INSTALL_DIR)lib
	@$(INSTALL) -c -m 644 $(EXT_LIB_DIR)* $(INSTALL_DIR)lib
	@$(INSTALL) -c -m 644 $(INC_DIR)*.h $(INSTALL_DIR)include/$(LIB_NAME)
	@$(RANLIB) $(INSTALL_DIR)lib/$(SLIB_FILE)

all: $(SLIB_FILE) tests tutorial

createFolders:
	@mkdir -p $(OBJ_DIR)

#-include $(DEPS)

# Dummy target to force rebuilds
FORCE:
