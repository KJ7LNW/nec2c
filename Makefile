#Makefile for nec2c     21 Aug 2003

SHELL = /bin/sh
CC    = gcc -Wall -O3 -g

objects = nec2c.o misc.o somnec.o

nec2c : $(objects)
	$(CC) -lm -lefence -o nec2c $(objects)

$(objects) : nec2c.h

nec2dx : 
	g77 -o nec2dx nec2dx.f

.PHONY : clean
clean  :
	-rm -f *.o *~

