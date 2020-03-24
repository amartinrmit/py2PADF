
VPATH = /usr/include:/scratch/amartin1/Work/codes/C-code

all : padfXcorr

padfXcorr : padfXcorr.o settings.o io.o interpolation.o qgrids.o fft1d.o correlation.o
	gcc -g -o padfXcorr padfXcorr.o settings.o io.o qgrids.o fft1d.o \
	correlation.o interpolation.o \
	-lm -L/usr/lib -lgsl -lgslcblas -lfftw3

padfXcorr.o  : padfXcorr.c
	gcc -g -c padfXcorr.c
 
settings.o  : settings.c
	gcc -g -c settings.c

io.o : io.c
	gcc -g -c io.c

qgrids.o : qgrids.c
	gcc -g -c qgrids.c

fft1d.o : fft1d.c
	gcc -g -c fft1d.c

interpolation.o : interpolation.c
	gcc -g -c interpolation.c

correlation.o : correlation.c
	gcc -g -c correlation.c

clean :
	rm padfXcorr padfXcorr.o settings.o io.o qgrids.o fft1d.o interpolation.o correlation.o
