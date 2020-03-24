/*
 * io.h      : part of padf
 * A.V.Martin March 2015
 *
 */



#ifndef IO_H
#define IO_H

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

void write_dbin( char * fname, double * array, int n);
void write_int_to_dbin( char * fname, int * array, int n);
void read_dbin( char * fname, double * array, int n);

long number_of_pixels_in_file( char * fname );
void allocate_and_read_binary_image( char * fname, long * n, double ** image);
void read_binary_image( char * fname, long n, double * image);


#endif
