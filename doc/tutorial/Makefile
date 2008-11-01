include config.mk

SRC = \
	genOscLFO1.cpp \
	genOscLFO2.cpp \
	genOscFMFeedback.cpp

TARGETS = ${SRC:.cpp=}

all: ${TARGETS}

%: %.cpp
	@echo CC $<
	@${CC} ${CFLAGS} -o $@ $@.cpp ${LDFLAGS}

clean:
	@rm -f *.o ${TARGETS}
