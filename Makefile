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

DEPFLAGS	=
ifneq ($(DEP_TRACK), 0)
	DEPFLAGS = -MMD -MF $(basename $@).dep
endif

#--------------------------------------------------------------------------
# Rules
#--------------------------------------------------------------------------

# Force these targets to always execute
.PHONY: clean cleanall external tests


# Build object file from C++ source
$(OBJ_DIR)%.o: %.cpp
	@echo CXX $< $@
	@$(CXX) -c $(CFLAGS) $< -o $@ $(DEPFLAGS)

-include $(wildcard $(OBJ_DIR)*.dep)


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
	@$(CXX) $(CFLAGS) -o $(BIN_DIR)$(*F) $@ $(SLIB_PATH) $(LDFLAGS)
ifneq ($(AUTORUN), 0)
	@cd $(BIN_DIR) && ./$(*F)
endif


# Remove active build configuration binary files
clean:
	$(call RemoveDir, $(OBJ_DIR))
	$(call RemoveDir, $(BIN_DIR))

# Remove all build configuration binary files
cleanall:
	@$(MAKE) clean BUILD_CONFIG=Release
	@$(MAKE) clean BUILD_CONFIG=Debug
	$(call RemoveDir, $(BUILD_DIR))


# Create file with settings for linking to external libraries
external:
	@printf "%b\n" "CPPFLAGS +=$(EXT_CPPFLAGS)\r\nLDFLAGS +=$(EXT_LDFLAGS)"> Makefile.external


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

test:
	@$(MAKE) tests/unitTests.cpp

# Archive repository
archive:
	$(eval $@_TMP := $(shell mktemp -d tmp.XXXXXXXXXX))
	@echo Creating archive, this may take some time...
	@echo Creating temporary export...
	@svn export --force . $($@_TMP)
	@echo Compressing...
	@cd $($@_TMP) && tar -czf ../$(LIB_NAME).tar.gz .
	@echo Compression complete.
	@$(RM) -R $($@_TMP)

createFolders:
	@mkdir -p $(OBJ_DIR)


# Dummy target to force rebuilds
FORCE:

