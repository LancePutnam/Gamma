# Customize to fit your system

# paths
PREFIX = /usr/local
CONFPREFIX = ${PREFIX}/etc
MANPREFIX = ${PREFIX}/share/man

# includes and libs
LIBS = -L${PREFIX}/lib -L/usr/lib


CFLAGS = -O3 -I. -I../include/ -I${PREFIX}/include -I/usr/include

# Linux/BSD
#LDFLAGS = ${LIBS} ../libgamma.a -lrt -lasound -ljack -lpthread -lportaudio -lsndfile

# Mac OS
LDFLAGS = ${LIBS} ../libgamma.a -lportaudio -lsndfile -framework AudioUnit -framework AudioToolbox -framework CoreAudio -framework Carbon

AR = ar cr
CC = g++
RANLIB = ranlib
