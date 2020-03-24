
#include "plotting.h"


void calculate_thplane( blarr_t * blarr, int nth, double theta, double rmax, int nl, double * output){

  int i, j, k;
  int nr = blarr->l[0].nq;
  double arg = theta;
  double sum;
  int alt;
  double leg;
  printf("DEBUG nr %d\n", nr);

  for(i=0;i<nr;i++){
    for(j=0;j<nr;j++){
  
      
      sum = 0.0;
  
      for(k=1;k<nl;k++){
	if ( (k%2)==1 ) {
	  alt = -1;
	} else {
	  alt = 1;
	}

	leg = gsl_sf_legendre_Pl ( k , cos(arg) );
	sum += blarr->l[k].data[i*nr+j] * leg * alt;
      } // k loop

      output[i*nr+j] = sum * 2 * PI; 
	//The correction for constant in the Legendre polynomial relation /// pow(2.0*PI,5.0);
    } //i loop
  } // j loop

}


void calculate_rplane( blarr_t * blarr, int nth, double r, double rmax, int nl, double * output){

  int i, j, k, i2;
  int nr = blarr->l[0].nq;
  double arg;
  double sum;
  int alt;
  double leg;


  i2 = (int) ((nr*r) / rmax);
  for(i=0;i<nr;i++){
    for(j=0;j<nth;j++){
      
      arg = j * 2 * PI / (double) nth;
      sum = 0.0;
      
      for(k=1;k<nl;k++){
	if ( (k%2)==1 ) {
	  alt = -1;
	} else {
	  alt = 1;
	}

	leg = gsl_sf_legendre_Pl ( k , cos(arg) );
	sum += blarr->l[k].data[i*nr+i2] * leg * alt;
      } // k loop

      output[i*nth+j] = sum * 2 * PI; 
	//The correction for constant in the Legendre polynomial relation /// pow(2.0*PI,5.0);
    } //i loop
  } // jloop

}

void calculate_r_eq_r_plane( blarr_t * blarr, int nth, double rmax, int nl, double * output){

  int i, j, k, i2;
  int nr = blarr->l[0].nq;
  double arg;
  double sum;
  int alt;
  double leg;


  for(i=0;i<nr;i++){
    for(j=0;j<nth;j++){
      
      arg = j * 2 * PI / (float) nth;
      sum = 0.0;
      
      for(k=1;k<nl;k++){
	if ( (k%2)==1 ) {
	  alt = -1;
	} else {
	  alt = 1;
	}

	leg = gsl_sf_legendre_Pl ( k , cos(arg) );
	sum += blarr->l[k].data[i*nr+i] * leg * alt;
      } // k loop

      output[i*nth+j] = sum * 2 * PI; 
	//The correction for constant in the Legendre polynomial relation /// pow(2.0*PI,5.0);
    } //i loop
  } // jloop

}


void calculate_rline( blarr_t * blarr, double theta, int nth, double r, double rmax, int nl, double * output){

  int i, j, k, i2;
  int nr = blarr->l[0].nq;
  int arg = theta;
  double sum;
  int alt;
  double leg;


  i2 = (int) ((nr*r) / rmax);
  for(i=0;i<nr;i++){
      
    sum = 0.0;
    
    for(k=1;k<nl;k++){
      if ( (k%2)==1 ) {
	alt = -1;
      } else {
	alt = 1;
      }
      
      leg = gsl_sf_legendre_Pl ( k , cos(arg) );
      sum += blarr->l[k].data[i*nr+i2] * leg * alt;
    } // k loop

    output[i] = sum * 2 * PI; 
	//The correction for constant in the Legendre polynomial relation /// pow(2.0*PI,5.0);
  } //i loop

}


void calculate_thline( blarr_t * blarr, int nth, double r, double r2, double rmax, int nl, double * output){

  int i, j, k, i2;
  int nr = blarr->l[0].nq;
  double arg;
  double sum;
  int alt;
  double leg;


  i = (int) ((nr*r) / rmax);
  i2 = (int) ((nr*r2) / rmax);
  printf("DEBUG i i2 %d %d %g %g\n", i, i2, r, r2);
  for(j=0;j<nth;j++){
    arg = j * 2 * PI / (float) nth;
    sum = 0.0;
    
    for(k=1;k<nl;k++){
      if ( (k%2)==1 ) {
	alt = -1;
      } else {
	alt = 1;
      }
      
      leg = gsl_sf_legendre_Pl ( k , cos(arg) );
      sum += blarr->l[k].data[i*nr+i2] * leg * alt;
    } // k loop

    output[j] = sum * 2 * PI; 
	//The correction for constant in the Legendre polynomial relation /// pow(2.0*PI,5.0);
  } //j loop

}
