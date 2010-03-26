#=========================================================================
# Gamma top level Makefile
#=========================================================================

include Makefile.config
SRCDIR = ./src

SRCS = 	$(SRCDIR)/Ambisonics.cpp\
	$(SRCDIR)/arr.cpp\
	$(SRCDIR)/AudioIO.cpp\
	$(SRCDIR)/Conversion.cpp\
	$(SRCDIR)/DFT.cpp\
	$(SRCDIR)/FFT_fftpack.cpp\
	$(SRCDIR)/fftpack++1.cpp\
	$(SRCDIR)/fftpack++2.cpp\
	$(SRCDIR)/File.cpp\
	$(SRCDIR)/Print.cpp\
	$(SRCDIR)/scl.cpp\
	$(SRCDIR)/Serialize.cpp\
	$(SRCDIR)/SoundFile.cpp\
	$(SRCDIR)/Sync.cpp\
	

OBJS = $(SRCS:.cpp=.o)

.cpp.o:
	@echo CC $<
	@$(CC) -c $(CFLAGS) -o $*.o $<

$(ALIB_FILE): $(OBJS)
	@echo AR $@
	@$(AR) $@ $(OBJS)
	@$(RANLIB) $@

$(SLIB_FILE): $(OBJS)
	@echo SH $@
	@$(CC) $(SLIBFLAGS) $(LFLAGS) -o $@ $(OBJS)

.PHONY: tests
tests: $(ALIB_FILE)
	@cd tests && make all

.PHONY: tutorial
tutorial: $(ALIB_FILE)
	@cd tutorial && make all

clean:
	@rm -f $(SRCDIR)/*.o *.$(ALIB_EXT) *.$(SLIB_EXT)
	@cd tests && make clean
	@cd tutorial && make clean

install: $(ALIB_FILE) $(SLIB_FILE)
	@$(INSTALL) -d $(PREFIX)/lib
	@$(INSTALL) -d $(PREFIX)/include/gamma
	$(INSTALL) -c -m 644 $(ALIB_FILE) $(PREFIX)/lib/$(ALIB_FILE)
	$(INSTALL) -c -m 644 $(SLIB_FILE) $(PREFIX)/lib/$(SLIB_FILE)
	$(INSTALL) -c -m 644 ./include/*.h $(PREFIX)/include/gamma
	$(RANLIB) $(PREFIX)/lib/$(ALIB_FILE)

#all: $(ALIB_FILE) $(SLIB_FILE) tests tutorial
all: $(ALIB_FILE) tests tutorial

