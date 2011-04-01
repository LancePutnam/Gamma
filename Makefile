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
	Recorder.cpp\
	SoundFile.cpp\
	Sync.cpp

OBJS = $(SRCS:.cpp=.o)
OBJS := $(addprefix $(OBJ_DIR), $(OBJS))
SRCS := $(addprefix $(SRC_DIR), $(SRCS))

CPPFLAGS	+= $(EXT_CPPFLAGS)
LDFLAGS		+= $(EXT_LDFLAGS)

CPPFLAGS	+= $(addprefix -I, $(INC_DIRS) $(RINC_DIRS))
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
	@$(CXX) -c $(CFLAGS) $< -o $@ -MMD -MF $(basename $@).dep

include $(wildcard $(OBJ_DIR)*.dep)


# Build static library
$(SLIB_PATH): createFolders $(OBJS)
	@echo AR $@
	@$(AR) $@ $(OBJS)
	@$(RANLIB) $@


# Build dynamic (shared) library
$(DLIB_FILE): createFolders external $(OBJS)
	@echo SH $@
	@$(CXX) $(DLIBFLAGS) $(LDFLAGS) -o $@ $(OBJS)


# Compile and run source files in examples/ and tests/ folders
EXEC_TARGETS = examples/%.cpp tests/%.cpp
.PRECIOUS: $(EXEC_TARGETS)
$(EXEC_TARGETS): $(SLIB_PATH) FORCE
	@$(CXX) $(CFLAGS) -o $(BIN_DIR)$(*F) $@ $(LDFLAGS) $(SLIB_PATH)
ifneq ($(AUTORUN), 0)
	@$(BIN_DIR)$(*F)
endif


# Remove active build configuration binary files
clean:
	@$(RM) $(OBJ_DIR)* $(OBJ_DIR) $(BIN_DIR)* $(BIN_DIR)


# Remove all build configuration binary files
cleanall:
	@$(MAKE) clean BUILD_CONFIG=release
	@$(MAKE) clean BUILD_CONFIG=debug
	@$(RM) $(BUILD_DIR)* $(BUILD_DIR)


# Create file with settings for linking to external libraries
external:
	@echo '\
CPPFLAGS += $(EXT_CPPFLAGS) \r\n\
LDFLAGS  += $(EXT_LDFLAGS) \
'> Makefile.external


# Clean and rebuild library
rebuild: clean $(SLIB_PATH)


# Install library into path specified by DESTDIR
# Include files are copied into DESTDIR/include/LIB_NAME and
# library files are copied to DESTDIR/lib
install: $(SLIB_PATH)
#	@echo 'INSTALL $(DESTDIR)'
	@$(INSTALL) -d $(DESTDIR)/lib
	@$(INSTALL) -d $(DESTDIR)/include/$(LIB_NAME)
	@$(INSTALL) -c -m 644 $(SLIB_PATH) $(DESTDIR)/lib
ifneq ($(EXT_LIB_COPY_DIR), )
	@$(INSTALL) -c -m 644 $(EXT_LIB_COPY_DIR)/* $(DESTDIR)/lib
endif
	@$(INSTALL) -c -m 644 $(INC_DIR)/*.h $(DESTDIR)/include/$(LIB_NAME)
	@$(RANLIB) $(DESTDIR)/lib/$(SLIB_FILE)


createFolders:
	@mkdir -p $(OBJ_DIR)


# Dummy target to force rebuilds
FORCE:

