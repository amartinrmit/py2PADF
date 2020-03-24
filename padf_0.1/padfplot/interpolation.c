
/*
 *  interpolation.c 
 *  A.V. Martin  March 2015
 *
 *  Use GSL to perform 2D interpolation; e.g. r-theta plots
 *
 */

#include "interpolation.h"



/*
 * 2D interpolation
 */ 
void interp_2D( interp_t * it ){

  int i,j,k;
  double out;
  gsl_interp ** gi = malloc( it->nx_in * sizeof(gsl_interp*) );
  gsl_interp_accel ** acc = malloc( it->nx_in * sizeof(gsl_interp_accel*));
  double ** din = malloc( it->nx_in * sizeof(double*) );
  double ** xin = malloc( it->nx_in * sizeof(double*) );
  double ** yin = malloc( it->nx_in * sizeof(double*) );

  // allocate the array of interpolation variables 
  for(i=0;i<it->nx_in;i++){
    gi[i] = gsl_interp_alloc ( gsl_interp_cspline, it->ny_in );
    acc[i] = gsl_interp_accel_alloc();
    din[i] = malloc( it->ny_in * sizeof(double) );
    xin[i] = malloc( it->ny_in * sizeof(double) );
    yin[i] = malloc( it->ny_in * sizeof(double) );
  }

  // loop over x values and create the splines
  
  for(i=0;i<it->nx_in;i++){
    for(j=0;j<it->ny_in;j++){
      din[i][j] = it->din[i*it->ny_in+j];
      xin[i][j] = it->xin[i*it->ny_in+j];
      yin[i][j] = it->yin[i*it->ny_in+j];
    }

    k = gsl_interp_init ( gi[i], yin[i], din[i], it->ny_in);

  }


  double * tempd  = malloc( it->nx_in * sizeof(double) );
  double * xclose = malloc( it->nx_in * sizeof(double) );
  gsl_interp * g = gsl_interp_alloc( gsl_interp_cspline, it->nx_in );
  gsl_interp_accel * accx = gsl_interp_accel_alloc();

  for(i=0;i<it->npix_out;i++){

    for(j=0;j<it->nx_in;j++){
      tempd[j] = gsl_interp_eval ( gi[j], yin[j], din[j], it->yout[i], acc[j]);
      x_search( it->nx_in, xin[j], it->xout[i], &(xclose[j]) );
    }
        

    k = gsl_interp_init ( g, xclose, tempd, it->ny_in );

    it->dout[i] = gsl_interp_eval ( g, xclose, tempd, it->xout[i], accx);
    
    k = gsl_interp_accel_reset( accx );
  }
  

  for(i=0;i<it->nx_in;i++){
    gsl_interp_free( gi[i] );
    gsl_interp_accel_free( acc[i] );
    free(din[i]);
  }
  free(din);
  free(gi);
  free(acc);
  gsl_interp_accel_free( accx );
  gsl_interp_free( g );
  free(tempd);
  free(xclose);
}


void x_search( int nx, double * xin, double xsample,  double * xclose ){

  int i;
  int istore=0;

  for (i=0;i<nx;i++){
    if(xsample < xin[i]){
      istore = i;
      break;
    }
  }
  
  *xclose = xin[istore];

}


void set_interp_t_input_array_sizes( interp_t * it, int nx_in, int ny_in ){

  it->nx_in = nx_in;
  it->ny_in = ny_in;
  it->npix_in = nx_in * ny_in;
}


void set_interp_t_output_array_sizes( interp_t * it, int nx_out, int ny_out ){

  it->nx_out = nx_out;
  it->ny_out = ny_out;
  it->npix_out = nx_out * ny_out;
}


void allocate_interp_t_arrays ( interp_t * it ){

  it->xin = malloc( it->npix_in * sizeof(double) );
  it->yin = malloc( it->npix_in * sizeof(double) );
  it->din = malloc( it->npix_in * sizeof(double) );

  it->xout = malloc( it->npix_out * sizeof(double) );
  it->yout = malloc( it->npix_out * sizeof(double) );
  it->dout = malloc( it->npix_out * sizeof(double) );

}


void free_interp_t_arrays ( interp_t * it ){

  free(it->xin);
  free(it->yin);
  free(it->din);

  free(it->xout);
  free(it->yout);
  free(it->dout);

}

void set_interp_t_input_xy_arrays( interp_t * it, double * xin, double * yin ){
  int i;
  
  for(i=0;i<it->npix_in;i++){
    it->xin[i] = xin[i];
    it->yin[i] = yin[i];
  }
}


void set_interp_t_input_data( interp_t * it, double * din ){
  int i;
  
  for(i=0;i<it->npix_in;i++){
    it->din[i] = din[i];
  }
}


void set_interp_t_output_arrays( interp_t * it, double * xout, double * yout ){
  int i;
  
  for(i=0;i<it->npix_out;i++){
    it->xout[i] = xout[i];
    it->yout[i] = yout[i];
  }

}

void xy_grids( int nx, int ny, double * x, double * y ){

  int i,j;
  for(i=0;i<nx;i++){
    for(j=0;j<ny;j++){
      x[i*ny+j] = i - nx/2;
      y[i*ny+j] = j - ny/2;
    }
  }

}

void xy_grids_scaled( int nx, int ny, double dmax, double * x, double * y ){

  int i,j;
  for(i=0;i<nx;i++){
    for(j=0;j<ny;j++){
      x[i*ny+j] = (i - nx/2) * dmax / (double)(nx/2);
      y[i*ny+j] = (j - ny/2) * dmax / (double)(nx/2);
    }
  }

}

// x y coordinate of each polar sampling point
void polar_xy_grids( int nr, int nth, double * rth_x, double * rth_y ){

  int i,j;

  for(i=0;i<nr;i++){
    for(j=0;j<nth;j++){
      rth_x[i*nth+j] = i * cos( 2 * PI * j / (double) nth );
      rth_y[i*nth+j] = i * sin( 2 * PI * j / (double) nth );
    }
  }

}
