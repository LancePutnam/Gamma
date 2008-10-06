# Customize to fit your system

# paths
PREFIX = /usr/local
CONFPREFIX = ${PREFIX}/etc
MANPREFIX = ${PREFIX}/share/man

# includes and libs
LIBS = -L${PREFIX}/lib -L/usr/lib

# Linux/BSD
CFLAGS = -g -I. -I./include/ -I${PREFIX}/include -I/usr/include
LDFLAGS = ${LIBS} -lsndfile -lportaudio -lfftwf -lm

AR = ar cr
CC = g++
RANLIB = ranlib
