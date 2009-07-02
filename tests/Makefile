include ../Makefile.config
LDFLAGS += ../libgamma.a
CFLAGS += -I../include/

SRC = \
	testAnalysis.cpp \
	testArray.cpp \
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
	testGenerator.cpp \
	testIpol.cpp \
	testMemory.cpp \
	testNode.cpp \
	testNoise.cpp \
	testOscillator.cpp \
	testRandom.cpp \
	testScalar.cpp \
	testTable.cpp \
	testTimer.cpp \
	testTypes.cpp \
	testUnits.cpp

TARGETS = ${SRC:.cpp=}

all: ${TARGETS}

%: %.cpp
	@echo CC $<
	@${CC} ${CFLAGS} -o $@ $@.cpp ${LDFLAGS}

clean:
	@rm -f *.o ${TARGETS}
