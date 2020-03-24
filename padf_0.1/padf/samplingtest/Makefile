VPATH = /usr/include:/home/amartin1/Work/codes/C-code

all : stest

stest : samplingtest.o 
	gcc -g -o stest samplingtest.o -lm -L/usr/lib -lgsl -lgslcblas -lfftw3

samplingtest.o  : samplingtest.c math.h
	gcc -g -c samplingtest.c

clean :
	rm stest samplingtest.o