#=========================================================================
# Gamma top level Makefile
#=========================================================================

include Makefile.config

SRCS = 	arr.cpp\
	Conversion.cpp\
	Domain.cpp\
	DFT.cpp\
	FFT_fftpack.cpp\
	fftpack++1.cpp\
	fftpack++2.cpp\
	Print.cpp\
	scl.cpp\
	Recorder.cpp\
	Scheduler.cpp\
	Timer.cpp

ifneq ($(NO_AUDIO_IO), 1)
	SRCS += AudioIO.cpp
endif

ifneq ($(NO_SOUNDFILE), 1)
	# needed for msys2 to compile sndfile.h
	ifeq ($(PLATFORM), windows)
		ifeq ($(MSYS_VERSION), 2)
			CPPFLAGS += -D __int64=int64_t
		endif
	endif
	SRCS += SoundFile.cpp
endif

#OBJS = $(SRCS:.cpp=.o)
#OBJS := $(addprefix $(OBJ_DIR), $(OBJS))
#SRCS := $(addprefix $(SRC_DIR), $(SRCS))

SRCS		:= $(addprefix $(SRC_DIR), $(SRCS))
OBJS		= $(addsuffix .o, $(basename $(notdir $(SRCS))))

CPPFLAGS	:= $(addprefix -I, $(INC_DIRS)) $(CPPFLAGS)


#--------------------------------------------------------------------------
# Rules
#--------------------------------------------------------------------------

gamma: $(LIB_PATH)

include Makefile.rules

# Force these targets to always execute
.PHONY: clean test


# Compile and run source files in RUN_DIRS
EXEC_TARGETS = $(addsuffix *.cpp, $(RUN_DIRS)) $(addsuffix *.c, $(RUN_DIRS)) $(addsuffix *.mm, $(RUN_DIRS))
.PRECIOUS: $(EXEC_TARGETS)
$(EXEC_TARGETS): $(LIB_PATH) FORCE
	$(CXX) $(ALL_CXXFLAGS) -o $(BIN_DIR)$(*F) $@ $(LIB_PATH) $(LDFLAGS)
ifneq ($(AUTORUN), 0)
	@cd $(BIN_DIR) && ./$(*F)
endif


# Remove active build configuration binary files
clean:
	$(call RemoveDir, $(OBJ_DIR))
	$(call RemoveDir, $(BIN_DIR))
	$(call RemoveDir, $(BUILD_DIR)lib/)
	$(call RemoveDir, $(BUILD_DIR)include/)
	$(call RemoveDir, $(BUILD_DIR))

# Clean and rebuild library
rebuild: clean $(LIB_PATH)


# Install library into path specified by DESTDIR
# Include files are copied into DESTDIR/include/LIB_NAME and
# library files are copied to DESTDIR/lib
install: $(LIB_PATH)
#	@echo 'INSTALL $(DESTDIR)'
	@$(INSTALL) -d $(DESTDIR)/lib
	@$(INSTALL) -d $(DESTDIR)/include/$(LIB_NAME)
	-@$(INSTALL) -m 644 $(LIB_PATH) $(DESTDIR)/lib
ifneq ($(EXT_LIB_COPY_DIR), )
	@$(INSTALL) -m 644 $(EXT_LIB_COPY_DIR)/* $(DESTDIR)/lib
endif
	@$(INSTALL) -m 644 $(INC_DIR)/*.h $(DESTDIR)/include/$(LIB_NAME)
#	@$(RANLIB) $(DESTDIR)/lib/$(LIB_FILE)

# Run unit tests
test:
	@$(MAKE) tests/unitTests.cpp RUN_DIRS=tests/

buildtest: test
	@for v in algorithmic analysis curves effects filter function io oscillator source spatial spectral synthesis synths techniques; do \
		$(MAKE) --no-print-directory examples/$$v/*.cpp RUN_DIRS=examples/$$v/ AUTORUN=0; \
	done


# Create/view API documentation
doc/html/index.html: doc/Doxyfile Gamma/*.h
	@if [ `which doxygen` ]; then \
		cd doc && doxygen Doxyfile && cd ..;\
	elif [ `which /Applications/Doxygen.app/Contents/Resources/doxygen` ]; then \
		cd doc && /Applications/Doxygen.app/Contents/Resources/doxygen Doxyfile && cd ..;\
	else \
		echo "Error: doxygen not found.";\
		echo "doxygen is required to create the documentation.";\
		printf "Please install it using ";\
		if [ `which apt-get` ]; then printf "\"sudo apt-get install doxygen\"";\
		elif [ `which port` ]; then printf "\"sudo port install doxygen\"";\
		elif [ `which brew` ]; then printf "\"brew install doxygen\"";\
		else printf "a package manager, e.g., apt-get (Linux), MacPorts or Homebrew (Mac OSX),";\
		fi;\
		printf " and try again. You can also create the documentation manually \
by downloading doxygen from www.doxygen.org and running it on the file $<.\n";\
		exit 127;\
	fi

docs: doc/html/index.html
ifeq ($(PLATFORM), linux)
	@xdg-open $< &
else ifeq ($(PLATFORM), macosx)
	@open $<
endif


