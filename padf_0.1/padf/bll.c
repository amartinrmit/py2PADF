/*
 *  bll.c
 *  A.V. Martin  March 2015
 *
 *
 */

#include "bll.h"


void allocate_correlation_t( correlation_t * corr, int nq, int nth ){

  corr->nq = nq;
  corr->nth = nth;
  corr->npix = nth*nq*nq;
  corr->data = malloc( nth*nq*nq*sizeof(double) );
}

void allocate_bl_t( bl_t * bl, int nq ){

  bl->nq = nq;
  bl->npix = nq*nq;
  bl->data = malloc( nq*nq*sizeof(double) );
  bl->noise = malloc( nq*nq*sizeof(double) );
  bl->filter = malloc( nq*sizeof(double) );

  int i;
  for (i=0;i<nq;i++){
    bl->filter[i] = 0.0;
  }
}

void allocate_blarr_t( blarr_t * blarr, int nl, int nq, int nlmin ){

  int i;
  blarr->nl = nl;
  blarr->nlmin = nlmin;
  blarr->l = malloc( nl*sizeof(bl_t) );

  for(i=0;i<nl;i++){
    allocate_bl_t( &(blarr->l[i]), nq );
  }
}

void allocate_blarr_t_sphBsamp( blarr_t * blarr, int nl, double rmax, double qmax ){

  int i;
  int nq;

  blarr->nl = nl;
  blarr->l = malloc( nl*sizeof(bl_t) );

  for(i=0;i<nl;i++){
    nq = sphB_samp_nmax( i, rmax, qmax );
    //printf("DEBUG blarr allocate nq %d\n", nq);
    allocate_bl_t( &(blarr->l[i]), nq );
  }
}


void allocate_fmat_t( fmat_t * fmat, int nl, int nth ){

  fmat->nl = nl/2;
  fmat->nth = nth;
  fmat->npix = nth*nl/2;
  fmat->data = malloc( nth*(nl/2)*sizeof(double) );
}

void allocate_noisedata_t( noisedata_t * nd, int nq ){
  nd->nq = nq;
  nd->n = nq*nq;
  nd->sigma = malloc( nq*nq* sizeof(double) );
}

void clean_up_correlation_t( correlation_t * corr ){
  free(corr->data);
}


void clean_up_bl_t( bl_t * bl ){
  free(bl->data);
  free(bl->noise);
  free(bl->filter);
}


void clean_up_blarr_t( blarr_t * blarr ){
  int i;
  for(i=0;i<blarr->nl;i++){
    clean_up_bl_t( &(blarr->l[i]) );
  }
}


void clean_up_fmat_t( fmat_t * fmat ){
  free(fmat->data);
}

void clean_up_noisedata_t( noisedata_t * nd ){
  free( nd->sigma);
}

void set_bl_t_blfilter( bl_t * bl, double blfilter, int blflag ){

  bl->blflag = blflag;
  bl->blfilter = blfilter;
}

void set_blarr_t_blfilter_uniform( blarr_t * bla, double blfilter, int blflag ){

  int i;
  for(i=0;i<bla->nl;i++){
    //bla->l[i].blflag = blflag;
    //bla->l[i].blfilter = blfilter;
    set_bl_t_blfilter( &(bla->l[i]), blfilter, blflag );    
  }
}


int sphB_samp_nmax( int l, double rmax, double qmax ){
  
  int i;
  int nmax = 1000000;
  int out = 0;
  double qln;
  double dl = l+0.5;

  for(i=1;i<nmax;i++){
    qln = gsl_sf_bessel_zero_Jnu ( dl, i );
    if ( qln > 2*PI*qmax*rmax ){
      out = i-2;
      break;
    }
  }

  if( out < 0 ){
    out = 0;
  }

  return out;

}



void fmat_calculate( fmat_t * fmat, double q, double q2, double  wl ){

  int i,j;
  double thetaq, thetaq2;
  double phi, arg;

  thetaq  = thetaq_calc( q,  wl );
  thetaq2 = thetaq_calc( q2, wl );
  

  for(i=0;i<fmat->nth;i++ ){
    for(j=0;j<fmat->nl;j++ ){
      
      phi = 2.0 * PI * i / (double) fmat->nth ;
      arg = cos(thetaq)*cos(thetaq2) + sin(thetaq)*sin(thetaq2)*cos(phi);
      if( arg > 1.0 ){
	arg = 1.0;
      }
      if( arg < -1.0 ){
	arg = -1.0;
      }

      fmat->data[i*fmat->nl+j] = (1.0/(4*PI))*gsl_sf_legendre_Pl ( j*2 , arg );
    }
  }

}



double thetaq_calc( double q, double wl ){

  double out = (PI/2.0) - asin( wl * q / (2.0) ); 
  return out;
}


void Blqq_calc( correlation_t * corr, blarr_t * bla, double wl, double rmax, double qmax ){

  int i,j,k,m;
  int ic, jc;
  int nq;
  int nth = corr->nth;
  double sum, qln1, qln2, dk;
  double thmax, th;
  fmat_t fmat;
  fmat_t inv;
  clock_t end, start, start2;

  double * x    = malloc( fmat.nl*sizeof(double) );
  double * temp = malloc( corr->nth*sizeof(double) );

  allocate_fmat_t( &fmat, bla->nl, corr->nth ); 
  allocate_fmat_t( &inv,  bla->nl, corr->nth ); 

  thmax = 2.0 * asin(qmax*wl/2.0);
  
  for(k=0;k<bla->nl;k++){
    
    nq = bla->l[k].nq;
    dk = k + 0.5;
    printf("DEBUG solving linear equations, l nq %d %d\n", k, nq );
    start = clock();
    for(i=0;i<bla->l[k].nq;i++){

      // if(k==2){
      //	printf("DEBUG q samples %g\n", gsl_sf_bessel_zero_Jnu ( dk, i+1 ) / rmax );
      //}

      for(j=0;j<bla->l[k].nq;j++){
	
	//printf("DEBUG k, i, j %d %d %d\n", k, i, j );

	qln1 =  gsl_sf_bessel_zero_Jnu ( dk, i+1 ) / (2*PI*rmax);
	qln2 =  gsl_sf_bessel_zero_Jnu ( dk, j+1 ) / (2*PI*rmax);

	fmat_calculate( &fmat, qln1, qln2, wl );

	// Only needed if svd is used
	//invert_fmat( &fmat, &inv );
    
	// nearest neighbour interpolation
	ic = round( (qln1/(qmax))*corr->nq );
	jc = round( (qln2/(qmax))*corr->nq );

	//th = 2.0 * asin(qln1*wl/2.0);
	//ic = round( (th/thmax)*corr->nq );
	//th = 2.0 * asin(qln2*wl/2.0);
	//jc = round( (th/thmax)*corr->nq );
	//printf("DEBUG ic jc %d %d %g %g\n", ic, jc, qln1, qln2);

	for(m=0;m<nth;m++){
	  temp[m] =  corr->data[ic*corr->nq*nth+jc*nth+m];
	}

	// QR version
	qr_solve( &fmat, temp, x );

	/* 	sum = 0.0;
	for(m=0;m<nth;m++){
	  // SVD version
	  // sum += inv.data[m*bla->nl+k] * corr->data[ic*corr->nq*nth+jc*nth+m];
	}
	*/

	//SVD version
	//bla->l[k].data[i*nq+j] = sum;

	//QR version
	if( (k%2)==0 ){
	  bla->l[k].data[i*nq+j] = x[k/2];
	} else {
	  bla->l[k].data[i*nq+j] = 0.0;
	}
      
      }
    }
    //printf("DEBUG l time %d %g %d\n", k, (double) (clock()-start)/CLOCKS_PER_SEC, bla->l[k].nq );
  }

  clean_up_fmat_t( &fmat );
  clean_up_fmat_t( &inv );
  free(x);
  free(temp);

}

/*
 *  constructs the basis functions using singular value decomposition
 */
void invert_fmat( fmat_t * fmat, fmat_t * inv ){

  int i,j,k;
  double *sout, *uout, *vout;
  gsl_matrix * A; 
  gsl_matrix * V; 
  gsl_vector * S;
  gsl_vector * work;

  // We assume that nth > nl

  int nth = fmat->nth, nl = fmat->nl;

  sout = malloc( nl*sizeof(double) );
  vout = malloc( nl*nl*sizeof(double) );
  uout = malloc( nl*nth*sizeof(double) );

  S = gsl_vector_alloc ( nl );
  work = gsl_vector_alloc ( nl );
   
  A = gsl_matrix_alloc ( nth, nl );
  V = gsl_matrix_alloc ( nl,  nl );

  for (i=0;i<nth;i++){
  for (j=0;j<nl;j++){

    gsl_matrix_set (A, i, j, fmat->data[i*nl+j]);

  }
  }

  int yes = gsl_linalg_SV_decomp_jacobi ( A, V, S);

  
  for (i=0;i<nl;i++){
    sout[i] = gsl_vector_get (S, i);
    //printf("%d %g;", i, sout[i] );
  }
  //printf("\n");
  //printf("DEBUG cond smax smin %g %g %g\n", sout[0]/sout[nl-1], sout[0], sout[nl-1]);
  

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
	if (sout[k] > sout[0]*1e-2 ){   //5e-2
	  sum2 += uout[i*n+k]*(1.0/sout[k])*vout[j*n+k];
	}
      }
      test[i*n+j] = sum;
      inv->data[i*n+j] = sum2;
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


/*
 *  constructs the basis functions using singular value decomposition
 */
void noise_matrix( int nth, int nl, double * mat, double * nmat ){

  int i,j,k;
  double *sout, *uout, *vout;
  gsl_matrix * A; 
  gsl_matrix * V; 
  gsl_vector * S;
  gsl_vector * work;

  // We assume that nth > nl

  //  int nth = fmat->nth, nl = nl;

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
    //printf("%d %g;", i, sout[i] );
  }
  //printf("\n");
  //printf("DEBUG cond smax smin %g %g %g\n", sout[0]/sout[nl-1], sout[0], sout[nl-1]);
  

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

  double sum = 0.0;
  // reconstruct the matrix and make the inverse
  for (i=0;i<m;i++){
    for (j=0;j<n;j++){
      sum = 0.0;
      for (k=0;k<n;k++){
	sum += uout[i*n+k]*sout[k]*vout[j*n+k] * uout[i*n+k]*sout[k]*vout[j*n+k];
      }
      nmat[i*n+j] = sum;
    }
  }


  free(sout);
  free(uout);
  free(vout);
  gsl_vector_free(work);
  gsl_vector_free(S);
  gsl_matrix_free(A);
  gsl_matrix_free(V);


}



/*
 *  constructs the basis functions using singular value decomposition
 */
void svd_solve( fmat_t * fmat, correlation_t * corr, noisedata_t * nd,
		int ic, int jc,
		int nq, double * output, double * noise_output ){

  int i,j,k;
  double *sout, *uout, *vout;
  gsl_matrix * A; 
  gsl_matrix * V; 
  gsl_vector * S;
  gsl_vector * work;

  // We assume that nth > nl

  int nth = fmat->nth, nl = fmat->nl;

  // a quick precondition
  /*
  for (i=0;i<nth;i++){
    if ((i<20)||(i>nth-20)){
      corr->data[ic*corr->nq*nth+jc*nth+i] *= 0.0;
      
      for (j=0;j<nl;j++){
	fmat->data[i*nl+j] *= 0.0;
      }
    }
  }
  */


  sout = malloc( nl*sizeof(double) );
  vout = malloc( nl*nl*sizeof(double) );
  uout = malloc( nl*nth*sizeof(double) );

  S = gsl_vector_alloc ( nl );
  work = gsl_vector_alloc ( nl );
   
  A = gsl_matrix_alloc ( nth, nl );
  V = gsl_matrix_alloc ( nl,  nl );

  for (i=0;i<nth;i++){
  for (j=0;j<nl;j++){

    gsl_matrix_set (A, i, j, fmat->data[i*nl+j]);

  }
  }

  int yes = gsl_linalg_SV_decomp_jacobi ( A, V, S);

  
  for (i=0;i<nl;i++){
    sout[i] = gsl_vector_get (S, i);
    //printf("%d %g;", i, sout[i] );
  }
  //printf("\n");
  //printf("DEBUG cond smax smin %g %g %g\n", sout[0]/sout[nl-1], sout[0], sout[nl-1]);
  

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
  double sum = 0.0;
  double sum2 = 0.0;
  
  double * udata = malloc( n * sizeof(double) );
  double * usig = malloc( n * sizeof(double) );

  for (j=0;j<n;j++){
    sum = 0.0;
    sum2 = 0.0;
    for (i=0;i<m;i++){
      sum += uout[i*nl+j] * corr->data[ic*corr->nq*nth+jc*nth+i];
      sum2 += uout[i*nl+j]*uout[i*nl+j]*nd->sigma[ic*corr->nq+jc]*nd->sigma[ic*corr->nq+jc];
    }
    udata[j] = sum;
    usig[j] = sum2;
    // if(j==0){
    //   printf("DEBUG svd_solve sum2 %g\n", sum2);
    // }
  }

  double filter;
  for (j=0;j<n;j++){
    if (sout[j] > sout[0]*5e-2 ){   //5e-2
      udata[j] *= 1.0 / sout[j];
      usig[j] *= 1.0 / (sout[j]*sout[j]);
    } else {
      udata[j] = 0.0;
      usig[j] = 0.0;
    }    
  }

  for (j=0;j<n;j++){
    sum = 0.0;
    sum2 = 0.0;
    for (i=0;i<n;i++){
      sum += vout[j*nl+i] * udata[i];
      sum2 += vout[j*nl+i] *vout[j*nl+i] * usig[i];
    }
    output[j] = sum;
   
    if (sum2>0.0){
      noise_output[j] = sqrt(sum2);
    }else{
      noise_output[j] = 0.0;
    }

  }


  double * test;
  double * test2;
  double * inv;
  test = malloc( m * n * sizeof(double) );
  inv = malloc( m * n * sizeof(double) );
  test2 = malloc( n * n * sizeof(double) );

  double* identity;
  identity = malloc( n * n * sizeof(double) );
  
 
  // reconstruct the matrix and make the inverse
  for (i=0;i<m;i++){
    for (j=0;j<n;j++){
      sum = 0.0;
      sum2 = 0.0;
      for (k=0;k<n;k++){
	sum += uout[i*n+k]*sout[k]*vout[j*n+k];
	if (sout[k] > sout[0]*1e-2 ){   //5e-2
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
  free(udata);
  free(usig);
  free(inv);


  free(sout);
  free(uout);
  free(vout);
  gsl_vector_free(work);
  gsl_vector_free(S);
  gsl_matrix_free(A);
  gsl_matrix_free(V);


}

/*
 *  constructs the basis functions using singular value decomposition
 */
void svd_solve_filter( fmat_t * fmat, correlation_t * corr, correlation_t * corrsig,
		       int ic, int jc,
		       int nq, double * output ){

  int i,j,k;
  double *sout, *uout, *vout;
  gsl_matrix * A; 
  gsl_matrix * V; 
  gsl_vector * S;
  gsl_vector * work;

  // We assume that nth > nl

  int nth = fmat->nth, nl = fmat->nl;

  // a quick precondition
  /*
  for (i=0;i<nth;i++){
    if ((i<20)||(i>nth-20)){
      corr->data[ic*corr->nq*nth+jc*nth+i] *= 0.0;
      
      for (j=0;j<nl;j++){
	fmat->data[i*nl+j] *= 0.0;
      }
    }
  }
  */


  sout = malloc( nl*sizeof(double) );
  vout = malloc( nl*nl*sizeof(double) );
  uout = malloc( nl*nth*sizeof(double) );

  S = gsl_vector_alloc ( nl );
  work = gsl_vector_alloc ( nl );
   
  A = gsl_matrix_alloc ( nth, nl );
  V = gsl_matrix_alloc ( nl,  nl );

  for (i=0;i<nth;i++){
  for (j=0;j<nl;j++){

    gsl_matrix_set (A, i, j, fmat->data[i*nl+j]);

  }
  }

  int yes = gsl_linalg_SV_decomp_jacobi ( A, V, S);

  
  for (i=0;i<nl;i++){
    sout[i] = gsl_vector_get (S, i);
    //printf("%d %g;", i, sout[i] );
  }
  //printf("\n");
  //printf("DEBUG cond smax smin %g %g %g\n", sout[0]/sout[nl-1], sout[0], sout[nl-1]);
  

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
  double sum = 0.0;
  double sum2 = 0.0;
  
  double * udata = malloc( n * sizeof(double) );
  double * usig = malloc( n * sizeof(double) );

  for (j=0;j<n;j++){
    sum = 0.0;
    sum2 = 0.0;
    for (i=0;i<m;i++){
      sum += uout[i*nl+j] * corr->data[ic*corr->nq*nth+jc*nth+i];
      sum2 += uout[i*nl+j] * corrsig->data[ic*corr->nq*nth+jc*nth+i]
	*uout[i*nl+j] * corrsig->data[ic*corr->nq*nth+jc*nth+i];
    }
    udata[j] = sum;
    usig[j]  = 0.001*sqrt(sum2);

  }

  double filter;
  for (j=0;j<n;j++){
    if (sout[j] > sout[0]*5e-2 ){   //5e-2
      filter = udata[j]*udata[j]/(udata[j]*udata[j] + usig[j]*usig[j]);
      udata[j] *= filter / sout[j];
    } else {
      udata[j] = 0.0;
    }
    // if(j==5){
    //  printf("DEBUG filter ic jc filter udata usig sout %d %d %g %g %g %g\n", ic, jc,
    //	     filter, udata[5], usig[5], sout[5]);
    //}
    
  }

  for (j=0;j<n;j++){
    sum = 0.0;
    for (i=0;i<n;i++){
      sum += vout[j*nl+i] * udata[i];
    }
    output[j] = sum;

  }


  double * test;
  double * test2;
  double * inv;
  test = malloc( m * n * sizeof(double) );
  inv = malloc( m * n * sizeof(double) );
  test2 = malloc( n * n * sizeof(double) );

  double* identity;
  identity = malloc( n * n * sizeof(double) );
  
 
  // reconstruct the matrix and make the inverse
  for (i=0;i<m;i++){
    for (j=0;j<n;j++){
      sum = 0.0;
      sum2 = 0.0;
      for (k=0;k<n;k++){
	sum += uout[i*n+k]*sout[k]*vout[j*n+k];
	if (sout[k] > sout[0]*1e-2 ){   //5e-2
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
  free(udata);
  free(usig);
  free(inv);


  free(sout);
  free(uout);
  free(vout);
  gsl_vector_free(work);
  gsl_vector_free(S);
  gsl_matrix_free(A);
  gsl_matrix_free(V);


}


/*
 *  constructs the basis functions using singular value decomposition
 */
void qr_solve( fmat_t * fmat, double * corr_data, double * output ){

  int i,j,k;
  double sum;
  int signum;
  gsl_matrix * A; 
  gsl_vector * x;
  gsl_vector * b;
  gsl_vector * tau;
  gsl_vector * res;
  gsl_vector * norm;
  gsl_permutation * p;
  gsl_matrix * q;
  gsl_matrix * r;
  gsl_matrix * r2;
  

  int yes;

  
  // We assume that nth > nl

  int nth = fmat->nth, nl = fmat->nl;

  A    = gsl_matrix_alloc ( nth, nl );
  x    = gsl_vector_alloc( nl );
  b    = gsl_vector_alloc( nl );
  tau  = gsl_vector_alloc( nl );
  norm = gsl_vector_alloc( nl );
  res  = gsl_vector_alloc( nth );
  p    = gsl_permutation_alloc( nl );

  q    = gsl_matrix_alloc ( nth, nth );
  r    = gsl_matrix_alloc ( nth, nl );
  r2   = gsl_matrix_alloc ( nl, nl );

  gsl_vector_set_zero( res) ;
  

  for (i=0;i<nth;i++){
  for (j=0;j<nl;j++){

    gsl_matrix_set (A, i, j, fmat->data[i*nl+j]);
   
  }
  }
  /*
  for (i=0;i<nth;i++){
    gsl_vector_set( b, i, corr_data[i]);
   
  }
  */

  yes = gsl_linalg_QRPT_decomp2( A, q, r, tau, p, &signum, norm );
  if( yes ){
    printf("GSL error: %s\n", gsl_strerror(yes) );
  }


  for (i=0;i<nl;i++){
    
    sum = 0.0;
    for(j=0;j<nth;j++){
      sum += gsl_matrix_get( q, j, i ) * corr_data[j];
    }

    gsl_vector_set( b, i, sum);  
  }

  for (i=0;i<nl;i++){
    for (j=0;j<nl;j++){
      gsl_matrix_set( r2, i, j, gsl_matrix_get(r, i, j ));
    }
  }


  //  yes = gsl_linalg_QR_lssolve( A, tau, b, x, res );
  yes = gsl_linalg_QRPT_Rsolve( r2, p, b, x );
  if( yes ){
    printf("GSL error: %s\n", gsl_strerror(yes) );
  }   


  for (j=0;j<nl;j++){
    output[j] = gsl_vector_get (x, j);
  }
  
  gsl_vector_free(x);
  gsl_vector_free(b);
  gsl_matrix_free(A);
  gsl_vector_free(tau);
  gsl_vector_free(res);
  gsl_vector_free(norm);
  gsl_permutation_free( p );
  gsl_matrix_free(q);
  gsl_matrix_free(r);
  gsl_matrix_free(r2);

}


/*
 *  constructs the basis functions using singular value decomposition
 */
void cgls_solve( fmat_t * fmat, double * corr_data, double * output, int rflag ){

  
  int i,j,k;
  int iter, niter;
  double sum, sum2;
  double alpha, beta;
  double r_error;

  double * xk;
  double * xk_p1;
  double * pk;
  double * pk_p1;
  double * rk;
  double * rk_p1;
  double * b;
  double * matrix;
  double * Ap;

  int nth = fmat->nth, nl = fmat->nl;

  xk     = malloc( nl   *sizeof(double) );
  xk_p1  = malloc( nl   *sizeof(double) );
  pk     = malloc( nl   *sizeof(double) );
  pk_p1  = malloc( nl   *sizeof(double) );
  rk     = malloc( nl   *sizeof(double) );
  rk_p1  = malloc( nl   *sizeof(double) );
  b      = malloc( nl   *sizeof(double) );
  matrix = malloc( nl*nl*sizeof(double) );
  Ap     = malloc( nl   *sizeof(double) );

  //set niter
  niter = 40;


  // create the square matrix
  for (i=0;i<nl;i++){
    for (j=0;j<nl;j++){

      sum = 0.0;
      for (k=0;k<nth;k++){
	sum += fmat->data[i*nl+k]*fmat->data[j*nl+k];
      }
      matrix[i*nl+j] = sum;
    }
  }

  // A_T * data
  for (j=0;j<nl;j++){
    sum = 0.0;
    for (k=0;k<nth;k++){
      sum += fmat->data[j*nl+k]*corr_data[k];
    }
    b[j] = sum;
  }

  // intialize cgls
  for (i=0;i<nl;i++){
    xk[i] = 0.0;
    rk[i] = b[i];
    pk[i] = rk[i];
  }


  for(iter=0;iter<niter;iter++){

    // Calculate Ap
    for (j=0;j<nl;j++){
      sum = 0.0;
      for (k=0;k<nl;k++){
	sum += matrix[j*nl+k]*pk[k];
      }
      Ap[j] = sum;
    } 

    // calcualte alpha
    sum = 0.0; sum2 = 0.0;
    for (i=0;i<nl;i++){
      sum  += rk[i]*rk[i];
      sum2 += pk[i]*Ap[i];
    }
    alpha = sum / sum2;

    // calculate xk_p1; rk_p1;
    for (i=0;i<nl;i++){
      xk_p1[i] = xk[i] + alpha * pk[i];
      rk_p1[i] = rk[i] - alpha * Ap[i];
    }

    // calcualte beta
    sum = 0.0; sum2 = 0.0;
    for (i=0;i<nl;i++){
      sum  += rk_p1[i]*rk_p1[i];
      sum2 += rk[i]*rk[i];
    }
    beta = sum / sum2;

    // calculate pk_p1
    for (i=0;i<nl;i++){
      pk_p1[i] = rk_p1[i] + beta * pk[i];
    }

    // reset xk, rk, pk
    for (i=0;i<nl;i++){
      xk[i] = xk_p1[i];
      rk[i] = rk_p1[i];
      pk[i] = pk_p1[i];
    }

    //residual error
    sum = 0.0; sum2 = 0.0;
    for (i=0;i<nl;i++){
      sum  += rk[i]*rk[i];
      sum2 += b[i]*b[i];
    }
    r_error = sum / sum2;
    
    if(rflag == 1){
      printf("Residual error: %d, %g\n", iter, r_error);
    }
  }


  for (i=0;i<nl;i++){
    output[i] = xk[i];
  }


  free(xk);
  free(xk_p1);
  free(pk);
  free(pk_p1);
  free(rk);
  free(rk_p1);
  free(b);
  free(matrix);
  free(Ap);

}



void bla_l_to_theta( blarr_t * bla, correlation_t * bth, correlation_t * bnoise,
		     int use_rl_filter, int noise_estimation_flag){
  
  int i,j,k,m;
  int nq = bth->nq;
  int nth = bth->nth;
  int nl = bla->nl;
  int nlmin = bla->nlmin;
  double arg;
  double sum, sum2;
  double alt;
  double leg;
  double rlfilter_tol = 0.05;
  double rlfactor=1.0;
  printf("DEBUG blar nq %d %d\n", bla->l[0].nq, bth->nq);
  printf("DEBUG l bla(r).filter[0] %g\n", bla->l[0].filter[125]);

  for(i=0;i<nq;i++){
    for(j=0;j<nq;j++){
  
      for(m=0;m<nth;m++){
	//arg = 2 * m / (double) nth - 1;
	arg = m * 2 * PI / (double) nth;
	sum = 0.0;

	for(k=nlmin;k<nl;k++){
	  if ( (k%2)==1 ) {
	    alt = -1;
	  } else {
	    alt = 1;
	  }

	  if (use_rl_filter == 1){
	    /*
	    printf("DEBUG i j k f0i %d %d %d %g \n", 
		   i, j, k, bla->l[0].filter[i]);
	    printf("DEBUG i j k f0j %d %d %d %g \n", 
		   i, j, k, bla->l[0].filter[j]);
	    printf("DEBUG i j k fki %d %d %d %g \n", 
		   i, j, k, bla->l[k].filter[i]);
	    printf("DEBUG i j k fkj %d %d %d %g \n", 
		   i, j, k, bla->l[k].filter[j]);
	    */
	    if ( (bla->l[0].filter[i]>0.0)
		 &&(bla->l[0].filter[j]>0.0)
		 &&(bla->l[k].filter[i] > bla->l[0].filter[i]*rlfilter_tol)
		 &&(bla->l[k].filter[j]> bla->l[0].filter[j]*rlfilter_tol)
		 ) {
	      rlfactor = 1.0 / (bla->l[k].filter[i]*bla->l[k].filter[j]);
	    } else {
	      rlfactor = 0.0;
	    }
	  }

	  leg = gsl_sf_legendre_Pl ( k , cos(arg) );
	  sum += bla->l[k].data[i*nq+j] * leg * alt * rlfactor;
	  //sum += bla->l[k].data[i*nq+j] * gsl_sf_legendre_Pl ( k , arg );
	}
	bth->data[i*nq*nth+j*nth+m] = sum * 2 * PI; //The correction for constant in the Legendre polynomial relation /// pow(2.0*PI,5.0);
	
	if (noise_estimation_flag ==1 ){
	  //	if (bnoise->data[0] != -1){
	  sum2 = 0.0;
	  for(k=nlmin;k<nl;k++){
	    if ( (k%2)==1 ) {
	      alt = -1;
	    } else {
	      alt = 1;
	    }
	    if (use_rl_filter == 1){
	      if ( (bla->l[0].filter[i]>0.0)
		   &&(bla->l[0].filter[j]>0.0)
		   &&(bla->l[k].filter[i] > bla->l[0].filter[i]*rlfilter_tol)
		   &&(bla->l[k].filter[j]> bla->l[0].filter[j]*rlfilter_tol)
		   ) {
		rlfactor = 1.0 / (bla->l[k].filter[i]*bla->l[k].filter[j]);
	      } else {
		rlfactor = 0.0;
	      }
	    }

	    leg = gsl_sf_legendre_Pl ( k , cos(arg) );
	    sum2 += bla->l[k].noise[i*nq+j] * leg * alt * bla->l[k].noise[i*nq+j] * leg * alt
	      *rlfactor*rlfactor;
	    // sum += bla->l[k].data[i*nq+j] * gsl_sf_legendre_Pl ( k , arg );
	  }
	  bnoise->data[i*nq*nth+j*nth+m] = sqrt(sum2) * 2 * PI;  
	 // }
	} //noise_estimation_flag - if

      } //m loop

    } //j loop
  } // i loop
  //printf("DEBUG bla_l_to_th i j m sum %d %d %d %g\n", i, j, m, sum);
}


void dsB_matrix( int nr, double rmax, int nmax, 
		 int l, double * mat){
  int i,j;
  double qln;
  double arg;
  double r;
  double factor;
  double jl1;
  double dl = l+0.5;

  for(i=0;i<nr;i++){
    r = i * rmax / (double) nr;

    for(j=0;j<nmax;j++){

      qln = gsl_sf_bessel_zero_Jnu ( dl, j+1 );
      arg = qln * r / rmax;
      jl1 = gsl_sf_bessel_jl ( l+1, qln );
      factor = sqrt(2 * PI) / ( pow(rmax,3) * jl1 * jl1 );

      mat[i*nmax+j] = gsl_sf_bessel_jl ( l, arg ) * factor;

    }
  }

}

void dsB_matrix_numerical_inv( int nr, double rmax, int nmax, 
			       int l, double * mat){

  double * tmpmat = malloc( nr*nmax*sizeof(double) );

  sphB_r_to_q_matrix( nr, rmax, nmax, l, tmpmat);
  // printf("DEBUG tmpmat : %g %g %g %g\n", tmpmat[(nr-3)*nmax+0], tmpmat[(nr-3)*nmax+1], tmpmat[(nr-3)*nmax+2], tmpmat[(nr-3)*nmax+3] );
  invert_qr_matrix( nmax, nr, tmpmat, mat, 0.02 );
  //printf("DEBUG mat : %g %g %g %g\n", mat[(nr-3)*nmax+0], mat[(nr-3)*nmax+1], mat[(nr-3)*nmax+2], mat[(nr-3)*nmax+3] );

  free(tmpmat);
}


void Bla_qr_transform( blarr_t * bla, blarr_t * bla_out, double qmax, double rmax, double * gfilter,
		       int gfilter_flag){

  int i,l;
  double * mat; 
  double * noisemat;
  double * rqmat;
  int nmax;
  int nq;
  int nr = bla_out->l[0].nq;

  for(l=0;l<bla->nl;l++){

    nq = bla->l[l].nq;
    mat = malloc( nq * nr * sizeof(double) );
    noisemat = malloc( nq * nr * sizeof(double) );
    rqmat = malloc( nq*nr*sizeof(double) );
 
    dsB_matrix( nr, rmax, nq, l, mat);
    //dsB_matrix_numerical_inv( nr, rmax, nq, l, mat);
    noise_matrix( nr, nq, mat, noisemat );
    sphB_r_to_q_matrix( nr, rmax, nq, l, rqmat);

    bl_qr_transform( nq, nr, mat, bla->l[l].data, bla_out->l[l].data );
    bl_qr_noise_transform( nq, nr, noisemat, bla->l[l].noise, bla_out->l[l].noise );
    calculate_bl_filter( nr, nq, rqmat, mat, bla_out->l[l].filter,
			 gfilter, gfilter_flag, rmax, qmax, l );
    if(l==0){ 
      printf("DEBUG l filter[0] %d %g\n", l, bla_out->l[0].filter[125]);
    }

    free(mat);
    free(noisemat);
    free(rqmat);
  }

}

//
//  spherical Bessel forward (r -> q) transform matrix
//
void sphB_r_to_q_matrix( int nr, double rmax, int nmax, 
		 int l, double * mat){

  int i,j;
  double qln;
  double arg;
  double r;
  double factor;
  double jl1;
  double dl = l + 0.5;

  for(i=0;i<nr;i++){
    r = (1.0*i) * rmax / (double) nr;

    for(j=0;j<nmax;j++){

      qln = gsl_sf_bessel_zero_Jnu ( dl, j+1 );
      //      arg = qln * r / rmax;
      arg = qln * i / (double) nr;
      jl1 = gsl_sf_bessel_jl ( l+1, qln );
      factor = sqrt(2.0 / PI); // .../ ( pow(rmax,3) * jl1 * jl1 );

      mat[i*nmax+j] = gsl_sf_bessel_jl ( l, arg ) * factor * r * r * rmax *1e30 /  (double) nr;
      if ( (i==nr-3)&&(j<4) ){
	//printf("DEBUG sphB : %d %g %d %g %g %g\n", i*nmax+j, gsl_sf_bessel_jl ( l, arg ), l, arg, qln, r/rmax );
      }
    }
  }
//  printf("DEBUG params : %g %g %g %g\n", factor, r, rmax, r*r );
//  printf("DEBUG tmpmat : %g %g %g %g\n", mat[(nr-3)*nmax+0], mat[(nr-3)*nmax+1], mat[(nr-3)*nmax+2], mat[(nr-3)*nmax+3] );
}

/*
 *  constructs the basis functions using singular value decomposition
 */
void invert_qr_matrix( int nmax, int nmax2, double * mat, double * inv, 
		       double thresh){

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
    //printf("DEBUG sout %d %g\n", i, sout[i]);
  }
  //printf("DEBUG cond smax smin %g %g %g\n", sout[0]/sout[nl-1], sout[0], sout[nl-1]);

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
	if (sout[k] > sout[0]*thresh){
	  sum2 += uout[i*n+k]*(1.0/sout[k])*vout[j*n+k];
	}
      }
      test[i*n+j] = sum;
      inv[i*n+j] = sum2;
    }
  }

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


void calculate_bl_filter(int nr, int nq, double * rqmat, double * mat, double * filter,
			 double * gfilter, int gfilter_flag, double rmax, double qmax, int l){

  int ir, ir2, iq;
  double sum, filtsum;
  double dl = l+0.5;
  double q, g;
  int ic;
 
  //  printf("DEBUG calcualte bl filter l nr %d %d\n", l, nr);
  for (ir=0;ir<nr;ir++){
    filtsum = 0.0;
    for(ir2=0;ir2<nr;ir2++){
      sum = 0.0;
      for(iq=0;iq<nq;iq++){

	if (gfilter_flag == 1){
	  q = gsl_sf_bessel_zero_Jnu ( dl, iq+1 ) / (2*PI*rmax);
	  //g = exp( - q * q / (gwid*gwid) );
	  ic = round( (q/(qmax))*nq );
	  g = gfilter[ic];
	} else {
	  g = 1.0;
	}
	//printf("DEBUG gwid q g %g %g %g\n", gwid, q, g);

	sum += rqmat[ir2*nq+iq]*mat[ir*nq+iq]*g;
      }
      filtsum += sum;
    }
    filter[ir] = filtsum;
    //   if(l==0){
    //   printf("DEBUG calcualte bl filter l ir filter[ir] %d %d %g\n", l, ir, filter[ir]);
    //  }
  }

}


void sphB_sampling_2( int l, int nq1, int nq2, double rmax, double * qx, double * qy ){

  int i,j;
  double qln1, qln2;
  double dl = l+0.5;

  for(i=0;i<nq1;i++){

    qln1 = gsl_sf_bessel_zero_Jnu ( dl, i+1 ) / rmax;

    for(j=0;j<nq2;j++){

      qln2 = gsl_sf_bessel_zero_Jnu ( dl, j+1 ) / rmax;
      
      qx[i*nq2+j] = qln1;
      qy[i*nq2+j] = qln2;

    }
  }

}


void bl_qr_transform( int nmax, int nr, double * mat, 
		      double * blin, double * blout ){

  int i,j;
  int iq, jq;
  double sum;

  for(i=0; i<nr; i++){
    for(j=0; j<nr; j++){
      
      sum = 0.0;
      for(iq=0; iq<nmax; iq++){
	for(jq=0; jq<nmax; jq++){
	  sum += blin[iq*nmax+jq] * mat[i*nmax+iq] * mat[j*nmax+jq];
	}
      }
      blout[i*nr+j] = sum / pow(2.0*PI,3.0);  // a correction for the normalization in the Discrete Spherical Bessel Transform
    }
  }
  
}

void bl_qr_noise_transform( int nmax, int nr, double * mat, 
			    double * blin, double * blout ){

  int i,j;
  int iq, jq;
  double sum;

  for(i=0; i<nr; i++){
    for(j=0; j<nr; j++){
      
      sum = 0.0;
      for(iq=0; iq<nmax; iq++){
	for(jq=0; jq<nmax; jq++){
	  sum += blin[iq*nmax+jq] *blin[iq*nmax+jq] * mat[i*nmax+iq] * mat[j*nmax+jq];
	}
      }
      if (sum>0.0){
	blout[i*nr+j] = sqrt(sum) / pow(2.0*PI,3.0);  // a correction for the normalization in the Discrete Spherical Bessel Transform
      } else {
	blout[i*nr+j] = 0.0;
      }
    }
  }
  
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
	//printf("DEBUG i j resamp %d %d %g\n", i, j, sum);
      }
    }
  }

}


/*
 *  constructs the basis functions using singular value decomposition
 */
void invert_resampling_matrix( int nmax, int nmax2, double * mat, double * inv,
			       double * noiseinv){

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
    //printf("DEBUG sout %d %g\n", i, sout[i]);
  }
  //printf("DEBUG cond smax smin %g %g %g\n", sout[0]/sout[nl-1], sout[0], sout[nl-1]);

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
  double sum3 = 0.0;
  

  // reconstruct the matrix and make the inverse
  for (i=0;i<m;i++){
    for (j=0;j<n;j++){
      sum = 0.0;
      sum2 = 0.0;
      sum3 = 0.0;
      for (k=0;k<n;k++){
	sum += uout[i*n+k]*sout[k]*vout[j*n+k];
	if (sout[k] > 0.5){
	  sum2 += uout[i*n+k]*(1.0/sout[k])*vout[j*n+k];
	  sum3 += (uout[i*n+k]*uout[i*n+k])*(1.0/(sout[k]*sout[k]))*(vout[j*n+k]*vout[j*n+k]);
	}
      }
      test[i*n+j] = sum;
      inv[i*n+j] = sum2;
      noiseinv[i*n+j] = sum3;
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

void Blqq_calc_fast( correlation_t * corr, correlation_t * corrsig, blarr_t * bla, 
		     double wl, double rmax, double qmax, int nfilterflag,
		     io_vars_t * iov, noisedata_t * nd  ){

  int i,j,k,m;
  int ic, jc;
  int nq;
  int nth = corr->nth;
  double sum, sum2, qln1, qln2, dk;
  fmat_t fmat;
  fmat_t inv;
  clock_t end, start, start2;
  double dbsum;
  double sigsum;

  resamp_t * rs;
  rs = malloc( bla->nl * sizeof(resamp_t) );
  resamp_t * rsinv;
  rsinv = malloc( bla->nl * sizeof(resamp_t) );
 
  blarr_t bltemp;
  allocate_blarr_t( &bltemp, bla->nl, bla->l[0].nq, bla->nlmin );
  for(i=0;i<bla->nl;i++){
    bltemp.l[i].blflag = bla->l[i].blflag;
    bltemp.l[i].blfilter = bla->l[i].blfilter;
  }

  printf("DEBUG Write blfilter test %g %g %d\n", bla->l[0].blfilter, bltemp.l[0].blfilter, bltemp.l[0].blflag);

  for(i=0;i<bla->nl;i+=2){
    rs[i].nmax1 = bla->l[i].nq;
    rs[i].nmax2 = bla->l[0].nq;
    rs[i].mat = malloc( bla->l[i].nq * bla->l[0].nq * sizeof(double) );
    start = clock();
    resampling_matrix( i, 0, rs[i].nmax1, rs[i].nmax2, qmax, rmax, rs[i].mat );
    printf("DEBUG calculated resampling matrix time %d %g\n", i, (double) (clock()-start)/CLOCKS_PER_SEC);
  }

  for(i=0;i<bla->nl;i+=2){
    rsinv[i].nmax1 = bla->l[i].nq;
    rsinv[i].nmax2 = bla->l[0].nq;
    rsinv[i].mat = malloc( bla->l[i].nq * bla->l[0].nq * sizeof(double) );
    rsinv[i].noisemat = malloc( bla->l[i].nq * bla->l[0].nq * sizeof(double) );
   
    start = clock();
    invert_resampling_matrix( rs[i].nmax1, rs[i].nmax2, rs[i].mat, rsinv[i].mat, rsinv[i].noisemat );
    //printf("DEBUG inverted resampling matrix time %d %g\n", i, (double) (clock()-start)/CLOCKS_PER_SEC);
  }

  // Need to define another bll array... with [0].nl

  allocate_fmat_t( &fmat, bla->nl, corr->nth ); 
  allocate_fmat_t( &inv,  bla->nl, corr->nth ); 
  
  double * x    = malloc( fmat.nl*sizeof(double) );
  double * temp = malloc( corr->nth*sizeof(double) );
  int rflag;


    
  nq = bla->l[0].nq;
  double * noise = malloc( nq*nq*sizeof(double) );
  dk = 0 + 0.5;
  //printf("DEBUG nq %d %d\n", k, nq );
  start = clock();
  for(i=0;i<bla->l[0].nq;i++){
    
    
    //printf("DEBUG q samples %g %g\n", gsl_sf_bessel_zero_Jnu ( dk, i+1 ) / (2*PI*rmax), qmax*i/(float)bla->l[0].nq );
    
    
    for(j=0;j<bla->l[0].nq;j++){
      
      //printf("DEBUG k, i, j %d %d %d\n", k, i, j );
      
      qln1 =  gsl_sf_bessel_zero_Jnu ( dk, i+1 ) / (2*PI*rmax);
      qln2 =  gsl_sf_bessel_zero_Jnu ( dk, j+1 ) / (2*PI*rmax);
      
      fmat_calculate( &fmat, qln1, qln2, wl );

      //SVD version
      invert_fmat( &fmat, &inv );

      // nearest neighbour interpolation
      ic = round( (qln1/(qmax))*corr->nq );
      jc = round( (qln2/(qmax))*corr->nq );
      //printf("DEBUG ic jc %d %d %g %g\n", ic, jc, qln1, qln2);
      
      if(nfilterflag==1){
	svd_solve_filter( &fmat, corr, corrsig, ic, jc, bla->l[0].nq, x );
      } else {
	svd_solve( &fmat, corr, nd, ic, jc, bla->l[0].nq, x, noise );
      }

      for(k=0;k<bla->nl;k+=2){
	bltemp.l[k].data[i*nq+j] = x[k/2];
	bltemp.l[k].noise[i*nq+j] = noise[k/2];
      }  

      for(k=1;k<bla->nl;k+=2){
	bltemp.l[k].data[i*nq+j] = 0.0;
	bltemp.l[k].noise[i*nq+j] = 0.0;
      }   


      // QR version
      /*
      for(m=0;m<nth;m++){
	temp[m] =  corr->data[ic*corr->nq*nth+jc*nth+m];
      }
      */

      //qr_solve( &fmat, temp, x );

      // CGLS Version
      /*
      if ((i==0)&&(j==0)){
	rflag = 1;
      } else {
	rflag = 0;
      }

      cgls_solve( &fmat, temp, x, rflag );
      */
      
      // SVD version
      /*
      for(k=0;k<bla->nl;k+=2){
	sum = 0.0;
	dbsum = 0.0;
	sigsum = 0.0;
	for(m=0;m<nth;m++){
	  sum += inv.data[m*fmat.nl+k/2] * corr->data[ic*corr->nq*nth+jc*nth+m];
	  sigsum += inv.data[m*fmat.nl+k/2] * corrsig->data[ic*corr->nq*nth+jc*nth+m]
	    *inv.data[m*fmat.nl+k/2] * corrsig->data[ic*corr->nq*nth+jc*nth+m];
	  dbsum += inv.data[m*fmat.nl+k/2];
	}
	bltemp.l[k].data[i*nq+j] = sum;
	noise[i*nq+j] = sqrt( sigsum );

	// noise filtering...
	if( (nfilterflag==1)&&( 3*noise[i*nq+j] > fabs(bltemp.l[k].data[i*nq+j]) ) ){
	  bltemp.l[k].data[i*nq+j] = 0.0;
	}
        */

	//printf("DEBUG corr/inv data sum %g %g\n", sum, dbsum);
      
      /*
       *  blfiltering
       */
      for(k=0;k<bla->nl;k+=2){
	if (bltemp.l[k].blflag==1){
	  //bltemp.l[k].data[i*nq+j] *=					\
	  //  exp(-(qln1-qln2)*(qln1-qln2)/(2.0*bltemp.l[k].blfilter*bltemp.l[k].blfilter));
	  bltemp.l[k].data[i*nq+j] *=					\
	    exp(-(qln1*qln1)/(4.0*bltemp.l[k].blfilter*bltemp.l[k].blfilter))
	    *exp(-(qln2*qln2)/(4.0*bltemp.l[k].blfilter*bltemp.l[k].blfilter));
	  //printf("DEBUG applied blfilter\n");
	}
      }
	/*  
	    printf("DEBUG ic jc delta_q filter %d %d %g %g %g %g %d %g\n", ic, jc, qln1, qln2, fabs(qln1-qln2), \
	    exp(-(qln1-qln2)*(qln1-qln2)/(2.0*bltemp.l[k].blfilter*bltemp.l[k].blfilter)),
	    bla->l[k].blflag, bla->l[k].blfilter ) ;
	*/
    
      
	/*
	  for(k=0;k<bla->nl;k+=2){
	  bltemp.l[k].data[i*nq+j] = x[k/2];
	  }   
	*/ 
   
    }
  }
  //printf("DEBUG l time %d %g %d\n", k, (double) (clock()-start)/CLOCKS_PER_SEC, bla->l[0].nq );
  
  char outname[1024];
  for(i=0;i<bla->nl;i+=2){
    set_fname_numbered( iov, "bltemp", ".dbin", i, outname );
    write_dbin( outname, bltemp.l[i].data, bla->l[0].nq*bla->l[0].nq );
  }


  // Add the resampled arrays into the bll terms.

  for(k=0;k<bla->nl;k++){
     nq = bla->l[k].nq;
     //printf("DEBUG k2 %d\n", (k%2) );
     if( (k%2)==0 ){
     
      for(i=0;i<bla->l[k].nq;i++){
	for(j=0;j<bla->l[k].nq;j++){
	  sum = 0.0;
	  sum2 = 0.0;
	  dbsum = 0.0;
	  for(ic=0;ic<bla->l[0].nq;ic++){
	    for(jc=0;jc<bla->l[0].nq;jc++){
	      sum += bltemp.l[k].data[ic*bla->l[0].nq+jc] * rsinv[k].mat[ ic*rsinv[k].nmax1+i] 
		* rsinv[k].mat[ jc*rsinv[k].nmax1+j ];
	      dbsum += fabs(bltemp.l[k].data[ic*bla->l[0].nq+jc]);

	      sum2 += bltemp.l[k].noise[ic*bla->l[0].nq+jc] *bltemp.l[k].noise[ic*bla->l[0].nq+jc]
		* rsinv[k].noisemat[ ic*rsinv[k].nmax1+i] * rsinv[k].noisemat[ jc*rsinv[k].nmax1+j ]
		* rsinv[k].noisemat[ ic*rsinv[k].nmax1+i] * rsinv[k].noisemat[ jc*rsinv[k].nmax1+j ];
	    }
	  }
	  bla->l[k].data[i*nq+j] = sum;
	  bla->l[k].noise[i*nq+j] = sqrt(sum2);
	  //if((j==2)&&(i==2)){printf("DEBUG sum %g %g\n", sum2, dbsum);}
	}
      }
     
     } else {
       for(i=0;i<bla->l[k].nq;i++){
	 for(j=0;j<bla->l[k].nq;j++){
	   bla->l[k].data[i*nq+j] = 0.0;
	   bla->l[k].noise[i*nq+j] = 0.0;
	 }
       }
     }

  }


  clean_up_fmat_t( &fmat );
  clean_up_fmat_t( &inv );
  free(temp);
  free(x);
  free(noise);
  for(i=0;i<bla->nl;i+=2){
    free(rs[i].mat);
    free(rsinv[i].mat);
    free(rsinv[i].noisemat);
  }
  free(rs);
  free(rsinv);
  clean_up_blarr_t( &bltemp );

}



void padf_calc( correlation_t * corr, correlation_t * corrsig, correlation_t * padf, 
		correlation_t * padf_noise,
		double rmax, double qmax,
		double wl, int nl, int nlmin, int blflag, double blfilter,
		int nfilterflag, io_vars_t * iov,
		noisedata_t * noisedata, int use_rl_filter, double * gfilter,
		int gfilter_flag, int noise_estimation_flag){

  int output_lmax = nl;
  char outname[1024];
  blarr_t blaq;
  blarr_t blar;
  clock_t start, end;

  allocate_blarr_t_sphBsamp( &blaq, nl, rmax, qmax );
  allocate_blarr_t( &blar, nl, padf->nq, nlmin );

  set_blarr_t_blfilter_uniform( &blaq, blfilter, blflag );
  printf("DEBUG Write blfilter test %g %d\n", blaq.l[0].blfilter, blaq.l[0].blflag);

  start = clock();
  // Blqq_calc( corr, &blaq, wl, rmax, qmax );
  Blqq_calc_fast( corr, corrsig, &blaq, wl, rmax, qmax, nfilterflag, iov,
		  noisedata );
  end = clock();
  printf("Blqq_calc took %ld seconds\n", (end-start)/CLOCKS_PER_SEC );

  
  int i;
  
  for(i=0;i<output_lmax;i+=2){
    set_fname_numbered( iov, "blqq", ".dbin", i, outname );
    if(i==0){printf("DEBUG writing first blqq\n");}
    write_dbin( outname, blaq.l[i].data, blaq.l[i].nq*blaq.l[i].nq );
    printf("DEBUG blqq i nq %d %d\n", i, blaq.l[i].nq);
  }
  

  //exit(0);

  start = clock();
  Bla_qr_transform( &blaq, &blar, qmax, rmax, gfilter, gfilter_flag );
  end = clock();
  printf("Bla_qr_transform took %ld seconds\n", (end-start)/CLOCKS_PER_SEC );

  
  for(i=0;i<output_lmax;i+=2){
    set_fname_numbered( iov, "blrr", ".dbin", i, outname );
    write_dbin( outname, blar.l[i].data, blar.l[i].nq*blar.l[i].nq );
    printf("DEBUG blrr i nq %d %d %g %d\n", i, blar.l[i].nq, blar.l[i].data[0], 
	   blar.l[i].nq*blar.l[i].nq);
  }
  
  printf("DEBUG l=0 blar.filter[0]%g\n", blar.l[0].filter[125]);
  

  start = clock();
  bla_l_to_theta( &blar, padf, padf_noise, use_rl_filter, noise_estimation_flag );
  end = clock();
  printf("bla_l_to_theta took %ld seconds\n", (end-start)/CLOCKS_PER_SEC );


  // strcpy( outname, "blar0.dbin" );
  //write_dbin( outname, blar.l[0].data, 128*128 );

  // output the l vs r filter
  strcpy( outname, iov->outpath );
  strcat( outname, "/" );
  strcat( outname, iov->tag );
  strcat( outname, "_r_vs_l_filter.txt" );
  FILE * f;
  f = fopen( outname, "w" );
  int ir, nr;
  nr = padf->nq;
  for(ir=0;ir<nr;ir++){
    //printf("%d %g ", ir, ir*rmax/(double) nr);
    fprintf(f, "%d %g ", ir, ir*rmax/(double) nr);
    for(i=0;i<nl;i+=2){
      //printf("%g ", blar[i]->filter[ir]);
      fprintf(f, "%g ", blar.l[i].filter[ir]);
    }
    //printf("\n");
    fprintf(f, "\n");
  }
  fclose(f);



  clean_up_blarr_t( &blar );
  clean_up_blarr_t( &blaq );
  printf("DEBUG Finished calculating padf\n");
}
