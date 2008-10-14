include config.mk

SRCDIR = ./src

SRC = \
	${SRCDIR}/Ambisonics.cpp  \
	${SRCDIR}/arr.cpp  \
	${SRCDIR}/AudioIO.cpp  \
	${SRCDIR}/DataFile.cpp  \
	${SRCDIR}/DFT.cpp  \
	${SRCDIR}/FFT_fftpack.cpp  \
	${SRCDIR}/fftpack.cpp  \
	${SRCDIR}/scl.cpp  \
	${SRCDIR}/SoundFile.cpp  \
	${SRCDIR}/Sync.cpp

OBJ = ${SRC:.cpp=.o}

.cpp.o:
	@echo CC $<
	@${CC} -c ${CFLAGS} -o $*.o $<

libgamma.a: ${OBJ}
	@echo AR $@
	@${AR} $@ ${OBJ}
	@${RANLIB} $@

.PHONY: tests

tests: libgamma.a
	@cd tests && make all

clean:
#	@rm -f ${SRCDIR}/*.o *.a
	@rm -f ${SRCDIR}/*.o
	@cd tests && make clean

all: libgamma.a tests
