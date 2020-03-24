/*
 *  correlation.h
 *  A.V. Martin  March 2015
 *
 *  Use fftw to perform autocorrelation of r-theta plot
 *
 */



#ifndef CORRELATION_H
#define CORRELATION_H


#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include "io.h"
#include "settings.h"
#include "fft1d.h"


void correlation_rtheta( settings * s, double * polar, double * corr );

void correlation_rtheta_faster( settings * s, double * polar, double * corr );

void correlation_X_rtheta( settings * s, double * polar, double * polar2, double * corr );

void correlation_X_rtheta_faster( settings * s, double * polar, double * polar2, double * corr );


#endif
