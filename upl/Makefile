#
# file -- Makefile 
# V-0.1
#
# Compile and Link for inversion
#

CC = /usr/bin/cc

SRC = fixed_threads_multi.c
#SRC = split_threads_multi.c
SRCH = globals.h

ENTRIES = fdfun.c \
          math_utils.c \
          kalman.c \
          text_utils.c \
          uncert.c
          
STAT_LIBS = -lm \
            -lpthread


#
# Definition of all program files.
#
PROGS = inversion

#
# TARGETS
#
all: ${PROGS}

install: all
	@echo nothing to install

clean:
	rm -f *.o *.gch lib/*.o lib/*.gch ${PROGS} ${GFILES} *.core

inversion: ${SRC} ${SRCH} ${ENTRIES}
	${CC} -g -o exec ${SRC} ${STAT_LIBS}


