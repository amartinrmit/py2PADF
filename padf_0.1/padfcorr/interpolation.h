/*
 *  interpolation.h
 *  A.V. Martin  March 2015
 *
 *  Use GSL to perform 2D interpolation; e.g. r-theta plots
 *
 */

#ifndef INTERPOLATION_H
#define INTERPOLATION_H



#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>
#include "qgrids.h"
#include "io.h"
#include "settings.h"

#define PI 4.0*atan(1.0)

typedef struct{

  int npix_in;
  int nx_in;
  int ny_in;
  double * xin;  // 2D array
  double * yin;  // 2D array
  double * din;

  int npix_out;
  int nx_out;
  int ny_out;
  double * xout;  // 2D array
  double * yout;  // 2D array
  double * dout;

} interp_t;


void interp_2D( interp_t * it );
void allocate_interp_t_arrays ( interp_t * it );
void free_interp_t_arrays ( interp_t * it );
void set_interp_t_input_array_sizes( interp_t * it, int nx_in, int ny_in );
void set_interp_t_output_array_sizes( interp_t * it, int nx_out, int ny_out );
void set_interp_t_input_xy_arrays( interp_t * it, double * xin, double * yin );
void set_interp_t_input_data( interp_t * it, double * din );
void set_interp_t_output_arrays( interp_t * it, double * xout, double * yout  );

// x_search function (return interpolated x-value)
void x_search( int nx, double * xin, double xsample,  double * xclose );

// grids with x and y coordinates of each pixel
// centred in middel of array.
//void xy_grids( int nx, int ny, double * x, double * y );
void xy_grids( int nx, int ny, double * x, double * y, int cx, int cy );
//void xy_grids_scaled( int nx, int ny, double dmax, double * x, double * y );
void xy_grids_scaled( int nx, int ny, double * x, double * y, int cx, int cy, double xscale, double yscale );
void polar_xy_grids( int nr, int nth, double * rth_x, double * rth_y );

// r-theta interpolation function
void polar_plot( settings * s, double * dp, double * polar_out );

void xy_grids_ewald( int nx, int ny, double * x, double * y,
		     double pw, double dz, double wl,
		     double qmax );

void polar_xy_grids_ewald( int nr, int nth, double * rth_x, double * rth_y, 
			   double pw, double dz, double wl,
			   double qmax );

double thetaq_calc( double qmax, double wl, int nr, int i );
#endif
