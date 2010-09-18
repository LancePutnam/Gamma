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
OBJS := $(addprefix $(OBJ_DIR)/, $(OBJS))

CFLAGS		+= $(addprefix -I, $(INC_DIRS) $(RINC_DIRS))
LFLAGS		:= $(addprefix -L, $(LIB_DIRS)) $(LFLAGS)
DLIB_FILE 	:= $(addprefix $(BIN_DIR)/, $(DLIB_FILE))
SLIB_FILE 	:= $(addprefix $(BIN_DIR)/, $(SLIB_FILE))

#--------------------------------------------------------------------------
# Targets
#--------------------------------------------------------------------------

# Build object file from C++ source
$(OBJ_DIR)/%.o: %.cpp
	@echo CC $< $@
	@$(CC) -c $(CFLAGS) $< -o $@

# Build static library
$(SLIB_FILE): createFolders $(OBJS)
	@echo AR $@
	@$(AR) $@ $(OBJS)
	@$(RANLIB) $@

# Build dynamic (shared) library
$(DLIB_FILE): createFolders $(OBJS)
	@echo SH $@
	@$(CC) $(DLIBFLAGS) $(LFLAGS) -o $@ $(OBJS)

.PHONY: tests
tests: $(SLIB_FILE)
	@make --directory $(TEST_DIR)

# Compile and run source files in examples/ folder
tests/%.cpp: $(SLIB_FILE) FORCE
	@$(CC) $(CFLAGS) -o $(BIN_DIR)/$(*F) $@ $(LFLAGS) $(SLIB_FILE)
	@./$(BIN_DIR)/$(*F)

.PHONY: tutorial
tutorial: $(SLIB_FILE)
	@make --directory $(TUT_DIR)

# Compile and run source files in examples/ folder
tutorial/%.cpp: $(SLIB_FILE) FORCE
	@$(CC) $(CFLAGS) -o $(BIN_DIR)/$(*F) $@ $(LFLAGS) $(SLIB_FILE)
	@$(BIN_DIR)/$(*F)

# Remove active build configuration binary files
.PHONY: clean
clean:
	@rm -rf $(BIN_DIR)/*
#	@find $(BIN_DIR) -type f ! -path '*.svn*' | xargs rm -f
#	@rm -f $(SRCDIR)/*.o *.$(SLIB_EXT) *.$(DLIB_EXT)
#	@cd tests && make clean
#	@cd tutorial && make clean

# Remove all build configuration binary files
.PHONY: cleanall
cleanall:
	@find $(BUILD_DIR) -type f ! -path '*.svn*' | xargs rm -f

install: $(SLIB_FILE) $(DLIB_FILE)
	@$(INSTALL) -d $(PREFIX)/lib
	@$(INSTALL) -d $(PREFIX)/include/gamma
	$(INSTALL) -c -m 644 $(SLIB_FILE) $(PREFIX)/lib/$(SLIB_FILE)
	$(INSTALL) -c -m 644 $(DLIB_FILE) $(PREFIX)/lib/$(DLIB_FILE)
	$(INSTALL) -c -m 644 ./include/*.h $(PREFIX)/include/gamma
	$(RANLIB) $(PREFIX)/lib/$(SLIB_FILE)

#all: $(SLIB_FILE) $(DLIB_FILE) tests tutorial
all: $(SLIB_FILE) tests tutorial

createFolders:
	@mkdir -p $(OBJ_DIR)

# Dummy target to force rebuilds
FORCE: