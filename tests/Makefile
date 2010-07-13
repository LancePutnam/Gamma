#=========================================================================
# Gamma tests Makefile
#=========================================================================

include ../Makefile.config
LFLAGS += -L../$(BIN_DIR)/ -l$(LIB_NAME)
CFLAGS += $(addprefix -I../, $(RINC_DIRS))

UTSRC = ut_arr.cpp \
	ut_fil.cpp \
	ut_gen.cpp \
	ut_ipl.cpp \
	ut_mem.cpp \
	ut_scl.cpp \
	utAccess.cpp \
	utContainers.cpp \
	utConversion.cpp \
	utFFT.cpp \
	utTypes.cpp

SRC = 	testAnalysis.cpp \
	testAudioChannels.cpp \
	testAudioSine.cpp \
	testContainers.cpp \
	testConversion.cpp \
	testDelay.cpp \
	testDFT.cpp \
	testEnvelope.cpp \
	testIpol.cpp \
	testNode.cpp \
	testNoise.cpp \
	testOscillator.cpp \
	testRandom.cpp \
	testRecurrenceMap.cpp \
	testScalar.cpp \
	testTable.cpp \
	testTimer.cpp \
	testTypes.cpp \
	testUnits.cpp

TARGETS = $(addprefix ../$(BIN_DIR)/, $(SRC:.cpp=) ut)
OBJS = $(addprefix ../$(OBJ_DIR)/, $(UTSRC:.cpp=.o))

#--------------------------------------------------------------------------
# Targets
#--------------------------------------------------------------------------
../$(OBJ_DIR)/%.o: %.cpp
	@echo CC $< $@
	@$(CC) -c $(CFLAGS) $< -o $@

all: $(TARGETS)

../$(BIN_DIR)/ut: $(OBJS)
	@echo CC $@
	@$(CC) $(CFLAGS) $^ -o $@ $(@F).cpp $(LFLAGS)

../$(BIN_DIR)/%: %.cpp
	@echo CC $<
	@$(CC) $(CFLAGS) -o $@ $^ $(LFLAGS)

