/*
 *  plotting.h
 *  A.V. Martin  March 2015
 *
 *
 */

#ifndef PLOTTING_H
#define PLOTTING_H



#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_fft_complex.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_linalg.h>
#include "io.h"
#include "settings.h"
#include "interpolation.h"
#include "bll.h"

void calculate_thplane( blarr_t * blarr, int nth, double theta, double rmax, int nl, double * output);
void calculate_rplane( blarr_t * blarr, int nth, double r, double rmax, int nl, double * output);
void calculate_r_eq_r_plane( blarr_t * blarr, int nth, double rmax, int nl, double * output);
void calculate_rline( blarr_t * blarr, double theta, int nth, double r, double rmax, int nl, double * output);
void calculate_thline( blarr_t * blarr, int nth, double r, double r2, double rmax,int nl,  double * output);

#endif
