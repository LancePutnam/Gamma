include Makefile.config
CFLAGS += -I./include/
SRCDIR = ./src

SRC = \
	${SRCDIR}/Ambisonics.cpp\
	${SRCDIR}/arr.cpp\
	${SRCDIR}/AudioIO.cpp\
	${SRCDIR}/Conversion.cpp\
	${SRCDIR}/DFT.cpp\
	${SRCDIR}/FFT_fftpack.cpp\
	${SRCDIR}/fftpack.cpp\
	${SRCDIR}/File.cpp\
	${SRCDIR}/scl.cpp\
	${SRCDIR}/SoundFile.cpp\
	${SRCDIR}/Sync.cpp\
	${SRCDIR}/Visual.cpp

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

.PHONY: tutorial
tutorial: libgamma.a
	@cd tutorial && make all

clean:
#	@rm -f ${SRCDIR}/*.o *.a
	@rm -f ${SRCDIR}/*.o
	@cd tests && make clean
	@cd tutorial && make clean

all: libgamma.a tests tutorial

