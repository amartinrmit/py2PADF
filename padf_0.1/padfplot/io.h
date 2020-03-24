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

typedef struct{

  char outpath[1024];
  char tag[1024];

} io_vars_t;



void write_dbin( char * fname, double * array, int n);
void write_int_to_dbin( char * fname, int * array, int n);
void read_dbin( char * fname, double * array, int n);

long number_of_pixels_in_file( char * fname );
void allocate_and_read_binary_image( char * fname, long * n, double ** image);
void read_binary_image( char * fname, long n, double * image);

void set_io_vars_t( io_vars_t * iovars, char * outpath, char * tag );
void set_fname( io_vars_t * iovars, char * label, char * fext, char * fname );
void set_fname_numbered( io_vars_t * iovars, char * label, char * fext, int i, char * fname  );

#endif
