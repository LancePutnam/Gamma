# Uncomment your system
#LINUX = 1
MAC = 1

# Customize to fit your system
PREFIX = /usr/local
CONFPREFIX = $(PREFIX)/etc
MANPREFIX = $(PREFIX)/share/man
LIBS = -L$(PREFIX)/lib -L/usr/lib
CFLAGS = -O3 -I. -I$(PREFIX)/include -I/usr/include
LDFLAGS = $(LIBS) -lsndfile -lportaudio -lm
AR = ar cr
CC = g++
RANLIB = ranlib

# OS dependent section
ifdef LINUX
	LDFLAGS += -lrt -lasound -ljack -lpthread
endif
ifdef MAC
	LDFLAGS += -framework AudioUnit -framework AudioToolbox -framework CoreAudio -framework Carbon
endif