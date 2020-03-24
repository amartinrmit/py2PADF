
/*
 *  samplingtest.c 
 *  A.V. Martin  April 2015
 *
 *  Test Bessel function resampling
 *
 */


#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_fft_complex.h>
#include <gsl/gsl_errno.h>

#define PI 4.0*atan(1.0)

void dsB_matrix_qtr( int nr, double qmax, double rmax, int nmax, 
		     int l, int l2, double * mat);
void dsB_matrix_rtq( int nq, double qmax, double rmax, int nr, 
		     int l, int l2, double * mat);
int sphB_samp_nmax( int l, double rmax, double qmax );
void resampling_matrix( int l, int l2, int nmax, int nmax2, 
			double qmax, double rmax, double * mat );
void invert_resampling_matrix( int nmax, int nmax2, double * fmat, double * inv );



void main(int argc, char * argv[]){

  int i,j,k;
  int nr, nr2;
  double sum;
  double rmax = 10.0;
  double qmax = 10.0;
  double qmax2;
  int l = 2;
  int l2 = 1;
  int nmax  = sphB_samp_nmax( l,  rmax, qmax );
  int nmax2 = sphB_samp_nmax( l2, rmax, qmax );
  double dl = l + 0.5;
  double dl2 = l2 + 0.5;
  nr = nmax;
  qmax =  nmax / rmax; // 2 * PI * nmax / rmax;
  qmax2 = nmax2 / rmax; // 2 * PI * nmax2 / rmax;

  double * func  = malloc( nr  * sizeof(double) );
  double r;
  double wid = nr/(4.0*2*PI*qmax);

  for(i=0;i<nr;i++){
    r = gsl_sf_bessel_zero_Jnu ( dl, i+1 ) / (2*PI*qmax);
    func[i] = exp( - r*r / (wid*wid));
    //printf("DEBUG i func %d %g\n", i, func[i]);
  }

  printf("DEBUG nmax etc %d %d %d %g\n", nmax, nmax2, nr, wid );


  double * forward = malloc( nr * nmax * sizeof(double) );
  double * forward2 = malloc( nr * nmax2 * sizeof(double) );
  double * back = malloc( nr * nmax2 * sizeof(double) );
  dsB_matrix_rtq( nmax, qmax, rmax, nr, l, l, forward );
  dsB_matrix_rtq( nmax2, qmax2, rmax, nr, l, l2, forward2 );
  dsB_matrix_qtr( nr, rmax, qmax, nmax, l, l, back );

  double * resamp = malloc( nmax * nmax2 * sizeof(double) );
  double * inv    = malloc( nmax * nmax2 * sizeof(double) );
  resampling_matrix( l, l2, nmax, nmax2, qmax, rmax, resamp );
  invert_resampling_matrix( nmax, nmax2, resamp, inv );

  // need to allocate new memory of output of maps
  // 1 for forward1, 1 for forward2, one for mapping of 1 to 2
  double * out1 = malloc( nmax * sizeof(double) );
  double * out2 = malloc( nmax2 * sizeof(double) );
  double * out3 = malloc( nmax * sizeof(double) );
  double * out4 = malloc( nmax2 * sizeof(double) );
  double * out5 = malloc( nr * sizeof(double) );


  for(i=0;i<nmax;i++){
    sum = 0.0;
    for(j=0;j<nr;j++){
      sum += forward[i*nr+j] * func[j] ;
    }
    out1[i] = sum;
  }

  for(i=0;i<nr;i++){
    sum = 0.0;
    for(j=0;j<nmax;j++){
      sum += back[i*nmax+j] * out1[j] ;
    }
    out5[i] = sum;
  }



  for(i=0;i<nmax2;i++){
    sum = 0.0;
    for(j=0;j<nr;j++){
      sum += forward2[i*nr+j] * func[j] ;
    }
    out2[i] = sum;
  }

   
  for(i=0;i<nmax;i++){
    sum = 0.0;
    for(j=0;j<nmax2;j++){
      sum += inv[j*nmax+i] * out2[j] ;
    }
    out3[i] = sum;
  }

  for(i=0;i<nmax2;i++){
    sum = 0.0;
    for(j=0;j<nmax;j++){
      sum += resamp[i*nmax+j] * out1[j] ;
    }
    out4[i] = sum;
  }
  
  
  int m,n;
  double q,q2,q3,q4;
  m = 1;
  n = 1;
  q = gsl_sf_bessel_zero_Jnu ( dl, m+1 );
  q2 = gsl_sf_bessel_zero_Jnu ( dl, n+1 );
  q4 = gsl_sf_bessel_zero_Jnu ( dl, nmax );
  sum = 0.0;
 
  for(j=0;j<nmax;j++){
    q3 = gsl_sf_bessel_zero_Jnu ( dl, j+1 );
    sum += 2.0 * PI * 2 * PI
      * (1.0/pow(gsl_sf_bessel_jl ( l+1, q ),2.0) )
      * (1.0/pow(gsl_sf_bessel_jl ( l+1, q3 ),2.0) )
      * gsl_sf_bessel_jl( l , q * q3 / q4   )
      * gsl_sf_bessel_jl( l, q2 * q3 / q4   )
      / pow( qmax*rmax, 3.0 );
  }
  printf("DEBUG m n sum %d %d %g\n", m, n, sum / (double) nmax);
  

  
  // check errors or print output or whatever
  printf("out1 out3\n");
  for(i=0;i<nmax;i++){
    printf("%g %g %g %g\n", gsl_sf_bessel_zero_Jnu ( dl, i+1 ) / rmax, out1[i], out3[i], out1[i]/out3[i]);
  }
  
  printf("\nout2 out4\n");
  for(i=0;i<nmax2;i++){
    printf("%g %g %g %g\n", gsl_sf_bessel_zero_Jnu ( dl2, i+1 ) / rmax, out2[i], out4[i], out2[i]/out4[i]);
  }
  
  printf("func out5\n");
  for(i=0;i<nr;i++){
    printf("%g %g %g %g\n", gsl_sf_bessel_zero_Jnu ( dl, i+1 ) / (2*PI*qmax), func[i], out5[i], out5[i]/func[i]);
  }
  
  printf("l l2 qmax qmax2 rmax nmax nmax2 %d %d %g %g %g %d %d\n", l, l2, qmax, qmax2, rmax, nmax, nmax2);
  printf("qmax*rmax %g %g\n", pow(qmax*rmax,1.0), 
	 sqrt(2*PI)/pow(gsl_sf_bessel_jl ( l+1,gsl_sf_bessel_zero_Jnu ( l, nmax )),2.0) );
  
  //free memory
  free(forward);
  free(forward2);
  free(func);
  free(resamp);
  free(inv);
  free(out1); free(out2); free(out3); free(out4);
}



// matrix that transforms q to r
void dsB_matrix_qtr( int nr, double rmax, double qmax, int nmax, 
		     int l, int l2, double * mat){
  int i,j;
  double qln;
  double arg;
  double r;
  double factor;
  double jl1;

  double dl = l + 0.5;
  double dl2 = l2 + 0.5;

  for(i=0;i<nr;i++){
    r = gsl_sf_bessel_zero_Jnu ( dl2, i+1 ) / (2 * PI * qmax);

    for(j=0;j<nmax;j++){

      qln = gsl_sf_bessel_zero_Jnu ( dl, j+1 );
      arg = r * (qln / rmax);
      jl1 = gsl_sf_bessel_jl ( l+1, qln );
      factor = sqrt(2 * PI) / ( pow(rmax,3) * jl1 * jl1 );
      //factor = sqrt(1.0) / ( pow(rmax,3) * jl1 * jl1 );

      mat[i*nmax+j] = gsl_sf_bessel_jl ( l, arg ) * factor ;

    }
  }

}

// matrix that transforms r to q
void dsB_matrix_rtq( int nq, double qmax, double rmax, int nr, 
		     int l, int l2, double * mat){
  int i,j;
  double rln;
  double arg;
  double q;
  double factor;
  double jl1;

  double dl = l + 0.5;
  double dl2 = l2 + 0.5;
  double qmax2p = qmax * 2 * PI;

  for(i=0;i<nq;i++){
    q = gsl_sf_bessel_zero_Jnu ( dl2, i+1 ) / rmax;

    for(j=0;j<nr;j++){

      rln = gsl_sf_bessel_zero_Jnu ( dl, j+1 );
      arg = q * (rln / qmax2p);
      jl1 = gsl_sf_bessel_jl ( l+1, rln );
      factor = sqrt(2 * PI) / ( pow(qmax2p,3) * jl1 * jl1 );
      //factor = sqrt(1.0) / ( pow(qmax,3) * jl1 * jl1 );

      mat[i*nr+j] = gsl_sf_bessel_jl ( l, arg ) * factor ;

    }
  }

}


int sphB_samp_nmax( int l, double rmax, double qmax ){
  
  int i;
  int nmax = 100000;
  int out = 0;
  double qln;
  double dl = l + 0.5;

  for(i=1;i<nmax;i++){
    qln = gsl_sf_bessel_zero_Jnu ( dl, i );
    //printf("DEBUG qln %g\n", qln);
    if ( qln > (2*PI*qmax*rmax) ){
      out = i-2;
      break;
    }
  }

  if( out < 0 ){
    out = 0;
  }
  printf("DEBUG out %d\n", out);

  return out;

}


void resampling_matrix( int l, int l2, int nmax, int nmax2, 
			double qmax, double rmax, double * mat ){

  int i,j,k;
  double sum;
  double q, q2;
  double r, factor, factor2;
  double jl1, jl2;
  double dl = l+0.5;
  double dl2 = l2+0.5;
  double qmax2p = qmax * 2 * PI;


  for(i=0;i<nmax2;i++){
    q2 = gsl_sf_bessel_zero_Jnu ( dl2, i+1 ) / rmax;
    
    for(j=0;j<nmax;j++){
      q = gsl_sf_bessel_zero_Jnu ( dl, j+1 );
      jl2 = gsl_sf_bessel_jl ( l+1, q );
      q *= 1.0/rmax;
      factor2 = sqrt(2 * PI) / ( pow(rmax,3) * jl2 * jl2 );

      sum = 0.0;
      for(k=0;k<nmax;k++){
	r = gsl_sf_bessel_zero_Jnu ( dl, k+1 );
	jl1 = gsl_sf_bessel_jl ( l+1, r );
	r *= 1.0/qmax2p;
	factor = sqrt(2 * PI) / ( pow(qmax2p,3) * jl1 * jl1 );

	sum += gsl_sf_bessel_jl ( l, q*r )*gsl_sf_bessel_jl ( l, q2*r )
	  * factor * factor2;
      }
      mat[i*nmax+j] = sum ;

      if (i==j){
	printf("DEBUG i j resamp %d %d %g\n", i, j, sum);
      }
    }
  }

}


/*
 *  constructs the basis functions using singular value decomposition
 */
void invert_resampling_matrix( int nmax, int nmax2, double * mat, double * inv ){

  int i,j,k;
  double *sout, *uout, *vout;
  gsl_matrix * A; 
  gsl_matrix * V; 
  gsl_vector * S;
  gsl_vector * work;

  // We assume that nth > nl

  int nth = nmax2, nl = nmax;

  sout = malloc( nl*sizeof(double) );
  vout = malloc( nl*nl*sizeof(double) );
  uout = malloc( nl*nth*sizeof(double) );

  S = gsl_vector_alloc ( nl );
  work = gsl_vector_alloc ( nl );
   
  A = gsl_matrix_alloc ( nth, nl );
  V = gsl_matrix_alloc ( nl,  nl );

  for (i=0;i<nth;i++){
  for (j=0;j<nl;j++){

    gsl_matrix_set (A, i, j, mat[i*nl+j]);

  }
  }

  int yes = gsl_linalg_SV_decomp_jacobi ( A, V, S);

  for (i=0;i<nl;i++){
    sout[i] = gsl_vector_get (S, i);
    printf("DEBUG sout %d %g\n", i, sout[i]);
  }
  printf("DEBUG cond smax smin %g %g %g\n", sout[0]/sout[nl-1], sout[0], sout[nl-1]);

  for (i=0;i<nth;i++){    
    for (j=0;j<nl;j++){
      uout[i*nl+j] = gsl_matrix_get (A, i, j) ;
    }
  }

  for (i=0;i<nl;i++){  
    for (j=0;j<nl;j++){
      vout[i*nl + j] = gsl_matrix_get (V, i, j);
    }
  }  
  
  int m = nth; 
  int n = nl;

  double * test;
  double * test2;
  test = malloc( m * n * sizeof(double) );
  test2 = malloc( n * n * sizeof(double) );

  double* identity;
  identity = malloc( n * n * sizeof(double) );
  
  double sum = 0.0;
  double sum2 = 0.0;

  // reconstruct the matrix and make the inverse
  for (i=0;i<m;i++){
    for (j=0;j<n;j++){
      sum = 0.0;
      sum2 = 0.0;
      for (k=0;k<n;k++){
	sum += uout[i*n+k]*sout[k]*vout[j*n+k];
	if (sout[k] > 0.5){
	  sum2 += uout[i*n+k]*(1.0/sout[k])*vout[j*n+k];
	}
      }
      test[i*n+j] = sum;
      inv[i*n+j] = sum2;
    }
  }

  /*  
  //make the identity
  for (i=0;i<n;i++){
    for (j=0;j<n;j++){
      identity[i*n+j] = 0.0;
      if (i==j){
	identity[i*n+j] = 1.0;
      }
     }
   }
   
  // test matrix
  sum = 0.0;
  for (i=0;i<m;i++){
    for (j=0;j<n;j++){
      sum += pow( fmat->data[i*n+j] - test[i*n+j], 2);
    }
  }
  //printf("matrix error %g\n", sum  );

  int ii,jj;
  for (i=0;i<n;i++){
    for (j=0;j<n;j++){
      sum = 0.0;
      for (k=0;k<m;k++){
	//sum += test[k*n+i]*inverse[k*n+j];
	sum += fmat->data[k*n+i]*inv->data[k*n+j];
	//sum += u[k*n+i]*u[k*n+test2];
      }
      test2[i*n+j] = sum;
    }
  }

  // test matrix
  //printf("(%g %g)\n ", test2[11*n+11] , identity[11*n+11]  );
  sum = 0.0;
  //printf("test2\n");
  for (i=0;i<n;i++){
    for (j=0;j<n;j++){
      sum += pow( identity[i*n+j] - test2[i*n+j], 2);
    }
  }
  printf("inverse error %g\n", sum  );
  */

  free(test);
  free(test2); 
  free(identity);


  free(sout);
  free(uout);
  free(vout);
  gsl_vector_free(work);
  gsl_vector_free(S);
  gsl_matrix_free(A);
  gsl_matrix_free(V);


}
