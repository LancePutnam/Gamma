include config.mk

SRC = \
	genEnvAD.cpp \
	genEnvDecay.cpp \
	genEnvSeg.cpp \
	genOscFMFeedback.cpp \
	genOscLFO1.cpp \
	genOscLFO2.cpp \
	genOscOsc.cpp \
	ioAudioDevice.cpp

TARGETS = ${SRC:.cpp=}

all: ${TARGETS}

%: %.cpp
	@echo CC $<
	@${CC} ${CFLAGS} -o $@ $@.cpp ${LDFLAGS}

clean:
	@rm -f *.o ${TARGETS}
