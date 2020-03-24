/*
 *  bll.h
 *  A.V. Martin  March 2015
 *
 *
 */

#ifndef BLL_H
#define BLL_H



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

typedef struct{ 
  int nq;
  int nth;
  int npix;
  double * data;
} correlation_t;

typedef struct{ 
  int nq;
  int npix;
  double * data;
  double * noise;
  double * filter;
  int blflag;
  double blfilter;
} bl_t;

typedef struct{ 
  int nl;
  int nlmin;
  bl_t * l;
} blarr_t;

typedef struct{ 
  int nth;
  int nl;
  int npix;
  double * data;
} fmat_t;

typedef struct{
  int nmax1;
  int nmax2;
  double * mat;
  double * noisemat;
} resamp_t;

typedef struct{
  int n;
  int nq;
  double * sigma; //standard deviation of noise on data
} noisedata_t;


// allocation routines
void allocate_correlation_t( correlation_t * corr, int nq, int nth );
void allocate_bl_t( bl_t * bl, int nq );
void allocate_blarr_t( blarr_t * blarr, int nl, int nq, int nlmin );
void allocate_blarr_t_sphBsamp( blarr_t * blarr, int nl, double rmax, double qmax );
void allocate_fmat_t( fmat_t * fmat, int nl, int nth );
void allocate_noisedata_t( noisedata_t * nd, int nq );

// clean up routines
void clean_up_correlation_t( correlation_t * corr );
void clean_up_bl_t( bl_t * bl );
void clean_up_blarr_t( blarr_t * blarr );
void clean_up_fmat_t( fmat_t * fmat );
void clean_up_noisedata_t( noisedata_t * nd );

void set_bl_t_blfilter( bl_t * bl, double blfilter, int blflag );
void set_blarr_t_blfilter_uniform( blarr_t * bla, double blfilter, int blflag );
int sphB_samp_nmax( int l, double rmax, double qmax );

// calculate Fmat
void fmat_calculate( fmat_t * fmat, double q, double q2, double  wl );
double thetaq_calc( double q, double wl );

// invert Fmat
void invert_fmat( fmat_t * fmat, fmat_t * inv );
void noise_matrix( int nth, int nl, double * mat, double * nmat );

void svd_solve( fmat_t * fmat, correlation_t * corr, noisedata_t * nd,
		int ic, int jc,
		int nq, double * output, double * noise_output );

void svd_solve_filter( fmat_t * fmat, correlation_t * corr, correlation_t * corrsig,
		       int ic, int jc,
		       int nq, double * output );

// QR solution
void qr_solve( fmat_t * fmat, double * corr_data, double * output );

void cgls_solve( fmat_t * fmat, double * corr_data, double * output, int rflag );

// calculate all Blqq matrices
void Blqq_calc( correlation_t * corr, blarr_t * bla, double wl, double rmax, double qmax );
void bla_l_to_theta( blarr_t * bla, correlation_t * bth, correlation_t * bnoise,
		     int use_rl_filter, int noise_estimation_flag);

void dsB_matrix( int nr, double rmax, int nmax, 
		 int l, double * mat);
void dsB_matrix_numerical_inv( int nr, double rmax, int nmax, 
			       int l, double * mat);
void sphB_r_to_q_matrix( int nr, double rmax, int nmax, 
			 int l, double * mat);

void Bla_qr_transform( blarr_t * bla, blarr_t * bla_out, double qmax, double rmax, 
		       double*gfilter, int gfilter_flag);

void sphB_sampling_2D( int l, int nq1, int nq2, double rmax, double * qx, double * qy );

void bl_qr_transform( int nmax, int nr, double * mat, 
		      double * blin, double * blout );
void bl_qr_noise_transform( int nmax, int nr, double * mat, 
			    double * blin, double * blout );

void invert_qr_matrix( int nmax, int nmax2, double * mat, double * inv, 
		       double thresh);
void calculate_bl_filter(int nr, int nq, double * rqmat, double * mat, double * filter,
			 double*gfilter, int gfilter_flag, double rmax, double qmax, int l);

void Blqq_calc_fast( correlation_t * corr, correlation_t * corrsig, blarr_t * bla, 
		     double wl, double rmax, double qmax, int nfilterflag,
		     io_vars_t * iov, noisedata_t * nd );

void padf_calc( correlation_t * corr, correlation_t * corrsig, correlation_t * padf, 
		correlation_t * padf_noise, double rmax, double qmax,
		double wl, int nl, int nlmin, int blflag, double blfilter, int nfilterflag, io_vars_t * iov,
		noisedata_t * nd, int use_rl_filter, double*gfilter, int gfilter_flag,
		int noise_estimation_flag);

void resampling_matrix( int l, int l2, int nmax, int nmax2, 
			double qmax, double rmax, double * mat );
void invert_resampling_matrix( int nmax, int nmax2, double * mat, double * inv, double * noiseinv );
#endif
