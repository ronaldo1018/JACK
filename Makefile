#CC=	g++
PROF=	-pg
#PROF=
#BASE=	-O3 -Wall
BASE=	-g  -Wall

CFLAGS= $(PROF) $(BASE)
LFLAGS= $(PROF) $(BASE) -lm -lpthread

HDR=	jack.h
SRC=	collect.c normal.c params.c proton.c queu.c table.c trace.c\
	delta.c tabmgr.c particle.c dump.c fluctuations.c fast.c  opencl.c cl_util.c
	
OBJ=	collect.o normal.o params.o proton.o queu.o table.o trace.o\
	delta.o tabmgr.o particle.o dump.o fluctuations.o fast.o  opencl.o cl_util.o

jack:	jack.o $(OBJ)
	$(CC) jack.o $(OBJ)  $(LFLAGS) -o jack -l OpenCL

search:	search.o $(OBJ)
	$(CC) search.o $(OBJ)  $(LFLAGS) -o search

hist:	hist.o
	$(CC) hist.o  $(LFLAGS) -o hist

convert: convert.o
	$(CC) convert.o  $(LFLAGS) -o convert

add:	add.o
	$(CC) add.o  $(LFLAGS) -o add

show:	show.o
	$(CC) show.o  $(LFLAGS) -o show

all:	jack search hist convert add show

test:	jack
	time ./jack -nmc 10000 -mean 1e5 -stdv 3e3

$(OBJ):	$(HDR)
jack.o:	jack.c $(HDR)
search.o: search.c $(HDR)

clean:
	rm -f $(OBJ) jack.o jack hist.o hist convert.o convert add.o add show.o show search.o search
