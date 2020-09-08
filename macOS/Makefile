#***************************************************************************   
#       user configuration
# set mpi = yes for the parallel version
# set debug = yes to debug 
DEBUG ?= no  
# end user configuration;
#***************************************************************************    

#SHELL = /bin/sh

CC = g++

ifeq ($(strip $(DEBUG)), yes)
    CFLAGS += -ggdb -DDEBUG
else
    CFLAGS += -O2
endif

CFLAGS += -Wall -m64 -std=c++11

#on mac, try ls /usr/lib/libz.*, there are some libz files in it.
#see https://stackoverflow.com/a/52511069 for installation of libz.a
STATICLIBS = libz.a

LIBS = -lm -lz

SRCSMAIN = main.cpp
SRCSFUNC = functions.cpp

HEADER = functions.h

.PHONY: all clean

all: nubeam-dedup

## note it's crutial to use -static-libgcc instead of static; 
nubeam-dedup: functions.o main.o
	$(CC) $(CFLAGS) functions.o main.o $(STATICLIBS) $(LIBS) -o nubeam-dedup

main.o: $(SRCSMAIN) $(HEADER) functions.o
	$(CC) $(CFLAGS) -c $(SRCSMAIN)

functions.o: $(SRCSFUNC) $(HEADER)
	$(CC) $(CFLAGS) -c $(SRCSFUNC)

clean: 
	@rm -v *.o
