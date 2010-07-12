#=========================================================================
# Gamma tutorials Makefile
#=========================================================================

include ../Makefile.config
LFLAGS += -L../$(BIN_DIR)/ -l$(LIB_NAME)
CFLAGS += -I../ -I../external/include

SRC =	effFlanger.cpp \
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
	ioSoundFile.cpp \
	vCrvBotanic.cpp \
	vCrvCircle.cpp \
	vCrvEpicycloid.cpp \
	vCrvHypocycloid.cpp \
	vCrvLissajous.cpp \
	vCrvSpiral.cpp

TARGETS = $(addprefix ../$(BIN_DIR)/, $(SRC:.cpp=))

#--------------------------------------------------------------------------
# Targets
#--------------------------------------------------------------------------
all: $(TARGETS)

../$(BIN_DIR)/%: %.cpp
	@echo CC $<
	@$(CC) $(CFLAGS) -o $@ $^ $(LFLAGS)
