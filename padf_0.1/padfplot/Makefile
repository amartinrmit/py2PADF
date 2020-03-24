
VPATH = /usr/include:/scratch/amartin1/Work/codes/C-code

all : padfplot

padfplot : padfplot.o settings.o io.o interpolation.o bll.o plotting.o
	gcc -g -o padfplot padfplot.o settings.o io.o \
	interpolation.o bll.o plotting.o \
	-lm -L/usr/lib -lgsl -lgslcblas -lfftw3

padfplot.o  : padfplot.c
	gcc -g -c padfplot.c

settings.o  : settings.c
	gcc -g -c settings.c

io.o : io.c
	gcc -g -c io.c

interpolation.o : interpolation.c
	gcc -g -c interpolation.c

bll.o : bll.c
	gcc -g -c bll.c

plotting.o : plotting.c
	gcc -g -c plotting.c

clean :
	rm padfplot padfplot.o settings.o io.o interpolation.o bll.o plotting.o
