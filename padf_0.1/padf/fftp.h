/*
 *  Routines for using the fftw library
 *  2D transforms
 *  Andrew V. Martin March 2015
 */


#ifndef FFTP_H
#define FFTP_H


#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <fftw3.h>
#include <time.h>

typedef struct{

  char fname[256];      // name of wisdom file
  int n;             // total pixels in 3D array   
  int nx;            // side length of array
  int ny;
  fftw_plan plan_forward;  // plan for the fft
  fftw_plan plan_back   ;  // plan for the fft
  fftw_complex * in;
  unsigned flags;

} fftvars;


void initialize_fftvars_2d( fftvars * v, int nx, int ny, char * fname );

void fftw_setup_2d( fftvars * v );

void fftw_clean_up( fftvars * v );

int fftw2d( fftvars * v, double * real, double * imag, 
	    double * real_out, double * imag_out, int n, const int sign );


#endif
