# Customize to fit your system

# paths
PREFIX = /usr/local
CONFPREFIX = ${PREFIX}/etc
MANPREFIX = ${PREFIX}/share/man

# includes and libs
LIBS = -L${PREFIX}/lib -L/usr/lib


CFLAGS = -g -I. -I../include/ -I${PREFIX}/include -I/usr/include

# Linux/BSD
#LDFLAGS = ${LIBS} ../libsynz.a -lportaudio -lsndfile -lfftw3f

# Mac OS
LDFLAGS = ${LIBS} ../libsynz.a -lportaudio -lsndfile -lfftw3f -framework AudioUnit -framework AudioToolbox -framework CoreAudio -framework Carbon

AR = ar cr
CC = g++
RANLIB = ranlib
