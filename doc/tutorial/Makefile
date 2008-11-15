include config.mk

SRC = \
	genEnvAD.cpp \
	genEnvDecay.cpp \
	genEnvSeg.cpp \
	genModAmp.cpp \
	genModFrqCM.cpp \
	genModFrqFeedback.cpp \
	genModFrqTransition.cpp \
	genNoise.cpp \
	genOscAddBeat.cpp \
	genOscAddPacket.cpp \
	genOscLFO1.cpp \
	genOscLFO2.cpp \
	genOscOsc.cpp \
	genOscSineD.cpp \
	ioAudioDevice.cpp \
	visLissajous.cpp \

TARGETS = ${SRC:.cpp=}

all: ${TARGETS}

%: %.cpp
	@echo CC $<
	@${CC} ${CFLAGS} -o $@ $@.cpp ${LDFLAGS}

clean:
	@rm -f *.o ${TARGETS}
