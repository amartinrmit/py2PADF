
VPATH = /usr/include:/scratch/amartin1/Work/codes/C-code

all : padf

padf : padf.o settings.o io.o interpolation.o bll.o
	gcc -g -o padf padf.o settings.o io.o \
	interpolation.o bll.o \
	-lm -L/usr/lib -lgsl -lgslcblas -lfftw3

#padf.o  : padf.c math.h
padf.o : padf.c
	gcc -g -c padf.c

settings.o  : settings.c
	gcc -g -c settings.c

io.o : io.c
	gcc -g -c io.c

interpolation.o : interpolation.c
	gcc -g -c interpolation.c

bll.o : bll.c
	gcc -g -c bll.c

clean :
	rm padf padf.o settings.o io.o interpolation.o bll.o
