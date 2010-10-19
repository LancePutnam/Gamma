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

CFLAGS		+= $(EXT_CFLAGS)
CFLAGS		+= $(addprefix -I, $(INC_DIRS) $(RINC_DIRS))
LFLAGS		+= $(EXT_LFLAGS)
LFLAGS		+= -L$(EXT_LIB_DIR)
DLIB_PATH 	:= $(addprefix $(BIN_DIR), $(DLIB_FILE))
SLIB_PATH 	:= $(addprefix $(BIN_DIR), $(SLIB_FILE))

DEPS		:= $(BUILD_DIR)depends

#--------------------------------------------------------------------------
# Targets
#--------------------------------------------------------------------------

# Force these targets to always execute
.PHONY: clean cleanall external tests

# Build object file from C++ source
$(OBJ_DIR)%.o: %.cpp
	@echo CC $< $@
	@$(CC) -c $(CFLAGS) $< -o $@

# Build static library
$(SLIB_PATH): createFolders external $(OBJS)
	@echo AR $@
	@$(AR) $@ $(OBJS)
	@$(RANLIB) $@

# Build dynamic (shared) library
$(DLIB_FILE): createFolders external $(OBJS)
	@echo SH $@
	@$(CC) $(DLIBFLAGS) $(LFLAGS) -o $@ $(OBJS)

$(DEPS): createFolders $(SRCS) FORCE
	@$(CC) -MM $(SRCS) -I./ > $(DEPS)

# Create file with settings for linking to external libraries
external:
	@echo '\
EXT_CFLAGS += $(EXT_CFLAGS) \r\n\
EXT_LFLAGS += $(EXT_LFLAGS) \
'> Makefile.external

# Compile and run source files in examples/ and tests/ folders
examples/%.cpp tests/%.cpp: $(SLIB_PATH) FORCE
	@$(CC) $(CFLAGS) -o $(BIN_DIR)$(*F) $@ $(LFLAGS) $(SLIB_PATH)
	@$(BIN_DIR)$(*F)

tests: $(SLIB_PATH)
	@$(MAKE) --directory $(TEST_DIR)

# Remove active build configuration binary files
clean:
	@rm -df $(OBJ_DIR)* $(BIN_DIR)*

# Remove all build configuration binary files
cleanall:
	@$(MAKE) clean BUILD_CONFIG=Release
	@$(MAKE) clean BUILD_CONFIG=Debug
	@rm -df $(BUILD_DIR)* $(BUILD_DIR)

# Clean and rebuild library
rebuild: clean $(SLIB_PATH)

install: $(SLIB_PATH)
	@$(INSTALL) -d $(INSTALL_DIR)
	@$(INSTALL) -d $(INSTALL_DIR)lib
	@$(INSTALL) -d $(INSTALL_DIR)include/$(INC_DIR)
	@$(INSTALL) -c -m 644 $(SLIB_PATH) $(INSTALL_DIR)lib
	@$(INSTALL) -c -m 644 $(EXT_LIB_DIR)* $(INSTALL_DIR)lib
	@$(INSTALL) -c -m 644 $(INC_DIR)*.h $(INSTALL_DIR)include/$(INC_DIR)
	@$(RANLIB) $(INSTALL_DIR)lib/$(SLIB_FILE)

all: $(SLIB_FILE) tests tutorial

createFolders:
	@mkdir -p $(OBJ_DIR)

#-include $(DEPS)

# Dummy target to force rebuilds
FORCE:
