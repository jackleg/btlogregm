# 
# Makefile for Tree-based Modeler
#
# Programmed by S.H.Lee
#
#

CC     = gcc
#CC=/opt/gcc-3.4.6-bc/bin/gcc
#CFLAGS = -g -Wall -fbounds-checking

#CFLAGS = -DDEBUG -Wall -Wunused -g
CFLAGS = -O3 -Wall -g

LOPT   = -lm
RFLAGS = -O3

BTREG  = btlogregdrvr

all : $(BTREG)

btlogregdrvr : btlogregdrvr.c btutil.c listmergesort.c listsort.c lnlist.c logregm.c Smemory.c getargs.c
	$(CC) $(CFLAGS) -o $@ btlogregdrvr.c btutil.c listmergesort.c listsort.c lnlist.c logregm.c Smemory.c getargs.c -lm

touch :
	touch *.c

clean : 
	rm -f *.o $(BTREG)
