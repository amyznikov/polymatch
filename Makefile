############################################################
#
# polymatch Makefile
# Generated by amyznikov Oct 18, 2016
#   from 'linux-gcc-executable' template
#
############################################################

SHELL = /bin/bash

TARGET = polymatch

all: $(TARGET)


cross   =
sysroot =
DESTDIR =
prefix  = /usr/local
bindir  = $(prefix)/bin
incdir  = $(prefix)/include
libdir  = $(prefix)/lib

INCLUDES+= -I.
SOURCES = $(wildcard *.c)
HEADERS = $(wildcard *.h)
MODULES = $(foreach s,$(SOURCES),$(addsuffix .o,$(basename $(s))))


# C preprocessor flags
CPPFLAGS=$(DEFINES) $(INCLUDES)

# C Compiler and flags
CC = $(cross)gcc
CFLAGS= -Wall -Wextra -Wno-unused-variable -Wno-unused-function -Wno-unused-parameter -O3 -g0

# Loader Flags And Libraries
LD=$(CC)
LDFLAGS = $(CFLAGS)

# STRIP = $(cross)strip --strip-all
STRIP = @echo "don't strip "

LDLIBS += -lm


#########################################


$(MODULES): $(HEADERS) Makefile
$(TARGET) : $(MODULES) Makefile
	$(LD) $(LDFLAGS)  $(MODULES) $(LDLIBS) -o $@

clean:
	$(RM) $(MODULES)

distclean: clean
	$(RM) $(TARGET)

install: $(TARGET) $(DESTDIR)/$(bindir)
	cp $(TARGET) $(DESTDIR)/$(bindir) && $(STRIP) $(DESTDIR)/$(bindir)/$(TARGET)

uninstall:
	$(RM) $(DESTDIR)/$(bindir)/$(TARGET)


$(DESTDIR)/$(bindir):
	mkdir -p $@
