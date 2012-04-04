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
	Scheduler.cpp\
	SoundFile.cpp\
	Sync.cpp

#OBJS = $(SRCS:.cpp=.o)
#OBJS := $(addprefix $(OBJ_DIR), $(OBJS))
#SRCS := $(addprefix $(SRC_DIR), $(SRCS))

SRCS		:= $(addprefix $(SRC_DIR), $(SRCS))
OBJS		= $(addsuffix .o, $(basename $(notdir $(SRCS))))

CPPFLAGS	+= $(addprefix -I, $(INC_DIRS) $(RINC_DIRS))
LDFLAGS		+= -L$(EXT_LIB_DIR)

#--------------------------------------------------------------------------
# Rules
#--------------------------------------------------------------------------

gamma: $(LIB_PATH)

include Makefile.rules

# Force these targets to always execute
.PHONY: clean cleanall external test


# Compile and run source files in examples/ and tests/ folders
EXEC_TARGETS = examples/%.cpp tests/%.cpp
.PRECIOUS: $(EXEC_TARGETS)
$(EXEC_TARGETS): $(LIB_PATH) FORCE
	@$(CXX) $(CFLAGS) -o $(BIN_DIR)$(*F) $@ $(SLIB_PATH) $(LDFLAGS)
#	@$(CXX) $(CFLAGS) -o $(BIN_DIR)$(*F) $@ -lGamma -Lbuild/lib $(LDFLAGS)
ifneq ($(AUTORUN), 0)
	@cd $(BIN_DIR) && ./$(*F)
endif


# Remove active build configuration binary files
clean:
	$(call RemoveDir, $(OBJ_DIR))
	$(call RemoveDir, $(BIN_DIR))


# Clean and rebuild library
rebuild: clean $(LIB_PATH)


# Install library into path specified by DESTDIR
# Include files are copied into DESTDIR/include/LIB_NAME and
# library files are copied to DESTDIR/lib
install: $(LIB_PATH)
#	@echo 'INSTALL $(DESTDIR)'
	@$(INSTALL) -d $(DESTDIR)/lib
	@$(INSTALL) -d $(DESTDIR)/include/$(LIB_NAME)
	@$(INSTALL) -c -m 644 $(LIB_PATH) $(DESTDIR)/lib
ifneq ($(EXT_LIB_COPY_DIR), )
	@$(INSTALL) -c -m 644 $(EXT_LIB_COPY_DIR)/* $(DESTDIR)/lib
endif
	@$(INSTALL) -c -m 644 $(INC_DIR)/*.h $(DESTDIR)/include/$(LIB_NAME)
#	@$(RANLIB) $(DESTDIR)/lib/$(LIB_FILE)

# Run unit tests
test:
	@$(MAKE) tests/unitTests.cpp

buildtest: test
	@for v in algorithmic curves effects filter function generator io spectral; do \
		$(MAKE) --no-print-directory examples/$$v/*.cpp AUTORUN=0; \
	done

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
