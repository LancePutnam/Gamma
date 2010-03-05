#=========================================================================
# Gamma tests Makefile
#=========================================================================

include ../Makefile.config
LFLAGS += ../libgamma.a
CFLAGS += -I../include/

UTSRC = ut_arr.cpp \
	ut_fil.cpp \
	ut_gen.cpp \
	ut_ipl.cpp \
	ut_mem.cpp \
	ut_scl.cpp \
	utAccess.cpp \
	utContainers.cpp \
	utConversion.cpp \
	utTypes.cpp

SRC = 	testAnalysis.cpp \
	testAudioBiquad.cpp \
	testAudioChannels.cpp \
	testAudioCombs.cpp \
	testAudioIO.cpp \
	testAudioLFO.cpp \
	testAudioMTOsc.cpp \
	testAudioNoise.cpp \
	testAudioPVoc.cpp \
	testAudioSampler.cpp \
	testAudioSine.cpp \
	testChaos.cpp \
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
	testScalar.cpp \
	testTable.cpp \
	testTimer.cpp \
	testTypes.cpp \
	testUnits.cpp

TARGETS = $(SRC:.cpp=) ut

.cpp.o:
	@echo CC $<
	@$(CC) -c $(CFLAGS) -o $*.o $<

all: $(TARGETS)

ut: $(UTSRC:.cpp=.o)
	@echo CC $@
	@$(CC) $(CFLAGS) $^ -o $@ $@.cpp $(LFLAGS)


%: %.cpp 
	@echo CC $<
	@$(CC) $(CFLAGS) -o $@ $@.cpp $(LFLAGS)

clean:
	@rm -f *.o $(TARGETS)
