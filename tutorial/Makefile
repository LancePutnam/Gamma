#=========================================================================
# Gamma tutorials Makefile
#=========================================================================

include ../Makefile.config
LFLAGS += ../libgamma.a
CFLAGS += -I../

SRC =	aEffFlanger.cpp \
	aIOFile.cpp \
	genEnvAD.cpp \
	genEnvDecay.cpp \
	genEnvSeg.cpp \
	genModAmp.cpp \
	genModFrqCM.cpp \
	genModFrqFeedback.cpp \
	genModFrqTransition.cpp \
	genNoise.cpp \
	genOscAddBeat.cpp \
	genOscAddBinBeat.cpp \
	genOscAddPacket.cpp \
	genOscLFO1.cpp \
	genOscLFO2.cpp \
	genOscOsc.cpp \
	genOscSineD.cpp \
	ioAudioDevice.cpp \
	vCrvBotanic.cpp \
	vCrvCircle.cpp \
	vCrvEpicycloid.cpp \
	vCrvHypocycloid.cpp \
	vCrvLissajous.cpp \
	vCrvSpiral.cpp

TARGETS = $(SRC:.cpp=)

all: $(TARGETS)

%: %.cpp
	@echo CC $<
	@$(CC) $(CFLAGS) -o $@ $@.cpp $(LFLAGS)

clean:
	@rm -f *.o $(TARGETS)
