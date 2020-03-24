
// Andrew V. Martin March 2015

#include "correlation.h"



void correlation_rtheta( settings * s, double * polar, double * corr ){

  int i, j, k;
  int nr = s->nr;
  int nth = s->nth;
  fftvars v;
  char fname[1024];
  double * in1d = malloc( s->nth * sizeof(double) );
  double * out1d = malloc( s->nth * sizeof(double) );
  double * imag = malloc( s->nth * sizeof(double) );

  double * in1d_b = malloc( s->nth * sizeof(double) );
  double * out1d_b = malloc( s->nth * sizeof(double) );
  double * imag_b = malloc( s->nth * sizeof(double) );

  double * tempr = malloc( s->nth * sizeof(double) );
  double * tempi = malloc( s->nth * sizeof(double) );

  for(j=0;j<s->nth;j++){
    imag[j] = 0.0;
    imag_b[j] = 0.0;
  }


  
  strcpy( fname, "wisdom.txt" );
  initialize_fftvars( &v, nth, fname );
  fftw_setup( &v );
  
  for(i=0;i<nr;i++){
    for(k=0;k<nr;k++){
    
      for(j=0;j<nth;j++){
	in1d[j] = polar[ i*nth+j];
	in1d_b[j] = polar[ k*nth+j];
	imag[j] = 0.0;
	imag_b[j] = 0.0;
      }
   
      fftw( &v, in1d, imag, out1d, imag, nth, 1 );
      fftw( &v, in1d_b, imag_b, out1d_b, imag_b, nth, 1 );
      
      for(j=0;j<nth;j++){
	tempr[j] = out1d[j] * out1d_b[j] + imag[j]*imag_b[j];
	tempi[j] = -out1d[j] * imag_b[j] + out1d_b[j]*imag[j];
      }

      fftw( &v, tempr, tempi, in1d, imag, nth, -1 );
      
      for(j=0;j<nth;j++){
	corr[ i*nth*nr+k*nth+j] = in1d[j]; //sqrt( in1d[j]*in1d[j] + imag[j]*imag[j]);
      }
      
    }
  }


  free(in1d);
  free(out1d);
  free(imag);
  free(in1d_b);
  free(out1d_b);
  free(imag_b);
  free(tempr);
  free(tempi);

}

void correlation_rtheta_faster( settings * s, double * polar, double * corr ){

  int i, j, k;
  int nr = s->nr;
  int nth = s->nth;
  fftvars v;
  char fname[1024];
  double * in1d = malloc( s->nth * sizeof(double) );
  double * out1d = malloc( s->nth * sizeof(double) );
  double * imag = malloc( s->nth * sizeof(double) );

  double * in1d_b = malloc( s->nth * sizeof(double) );
  double * out1d_b = malloc( s->nth * sizeof(double) );
  double * imag_b = malloc( s->nth * sizeof(double) );

  double * tempr = malloc( s->nth * sizeof(double) );
  double * tempi = malloc( s->nth * sizeof(double) );

  double * fpolar_real = malloc( s->nth * s->nr * sizeof(double) );
  double * fpolar_imag = malloc( s->nth * s->nr * sizeof(double) );


  for(j=0;j<s->nth;j++){
    imag[j] = 0.0;
    imag_b[j] = 0.0;
  }

  for(i=0;i<nr;i++){ 
    for(j=0;j<nth;j++){
      fpolar_real[i*nth+j] = 0.0;
      fpolar_imag[i*nth+j] = 0.0;
    }
  }


  strcpy( fname, "wisdom.txt" );
  initialize_fftvars( &v, nth, fname );
  fftw_setup( &v );


  for(i=0;i<nr;i++){ 
    for(j=0;j<nth;j++){
      in1d[j] = polar[ i*nth+j];
      imag[j] = 0.0;
    }
  
    fftw( &v, in1d, imag, out1d, imag, nth, 1 );

    for(j=0;j<nth;j++){
      fpolar_real[ i*nth+j] = out1d[j];
      fpolar_imag[ i*nth+j] = imag[j];
    }

  }

  
  
  for(i=0;i<nr;i++){
    for(k=0;k<nr;k++){
    
      for(j=0;j<nth;j++){
	out1d[j] = fpolar_real[ i*nth+j];
	out1d_b[j] = fpolar_real[ k*nth+j];
	imag[j] = fpolar_imag[ i*nth+j];
	imag_b[j] = fpolar_imag[ k*nth+j];
      }
   
      for(j=0;j<nth;j++){
	tempr[j] = out1d[j] * out1d_b[j] + imag[j]*imag_b[j];
	tempi[j] = -out1d[j] * imag_b[j] + out1d_b[j]*imag[j];
      }

      fftw( &v, tempr, tempi, in1d, imag, nth, -1 );
      
      for(j=0;j<nth;j++){
	corr[ i*nth*nr+k*nth+j] = in1d[j]; //sqrt( in1d[j]*in1d[j] + imag[j]*imag[j]);
      }
      
    }
  }


  free(fpolar_real);
  free(fpolar_imag);
  free(in1d);
  free(out1d);
  free(imag);
  free(in1d_b);
  free(out1d_b);
  free(imag_b);
  free(tempr);
  free(tempi);

}



void correlation_X_rtheta( settings * s, double * polar, double * polar2, double * corr ){

  int i, j, k;
  int nr = s->nr;
  int nth = s->nth;
  fftvars v;
  char fname[1024];
  double * in1d = malloc( s->nth * sizeof(double) );
  double * out1d = malloc( s->nth * sizeof(double) );
  double * imag = malloc( s->nth * sizeof(double) );

  double * in1d_b = malloc( s->nth * sizeof(double) );
  double * out1d_b = malloc( s->nth * sizeof(double) );
  double * imag_b = malloc( s->nth * sizeof(double) );

  double * tempr = malloc( s->nth * sizeof(double) );
  double * tempi = malloc( s->nth * sizeof(double) );

  for(j=0;j<s->nth;j++){
    imag[j] = 0.0;
    imag_b[j] = 0.0;
  }


  
  strcpy( fname, "wisdom.txt" );
  initialize_fftvars( &v, nth, fname );
  fftw_setup( &v );
  
  for(i=0;i<nr;i++){
    for(k=0;k<nr;k++){
    
      for(j=0;j<nth;j++){
	in1d[j] = polar[ i*nth+j];
	in1d_b[j] = polar2[ k*nth+j];
	imag[j] = 0.0;
	imag_b[j] = 0.0;
      }
   
      fftw( &v, in1d, imag, out1d, imag, nth, 1 );
      fftw( &v, in1d_b, imag_b, out1d_b, imag_b, nth, 1 );
      
      for(j=0;j<nth;j++){
	tempr[j] = out1d[j] * out1d_b[j] + imag[j]*imag_b[j];
	tempi[j] = -out1d[j] * imag_b[j] + out1d_b[j]*imag[j];
      }

      fftw( &v, tempr, tempi, in1d, imag, nth, -1 );
      
      for(j=0;j<nth;j++){
	corr[ i*nth*nr+k*nth+j] = in1d[j]; //sqrt( in1d[j]*in1d[j] + imag[j]*imag[j]);
      }
      
    }
  }


  free(in1d);
  free(out1d);
  free(imag);
  free(in1d_b);
  free(out1d_b);
  free(imag_b);
  free(tempr);
  free(tempi);

}




void correlation_X_rtheta_faster( settings * s, double * polar, double * polar2, double * corr ){

  int i, j, k;
  int nr = s->nr;
  int nth = s->nth;
  fftvars v;
  char fname[1024];
  double * in1d = malloc( s->nth * sizeof(double) );
  double * out1d = malloc( s->nth * sizeof(double) );
  double * imag = malloc( s->nth * sizeof(double) );

  double * in1d_b = malloc( s->nth * sizeof(double) );
  double * out1d_b = malloc( s->nth * sizeof(double) );
  double * imag_b = malloc( s->nth * sizeof(double) );

  double * tempr = malloc( s->nth * sizeof(double) );
  double * tempi = malloc( s->nth * sizeof(double) );

  double * fpolar_real = malloc( s->nth * s->nr * sizeof(double) );
  double * fpolar_imag = malloc( s->nth * s->nr * sizeof(double) );
  double * fpolar2_real = malloc( s->nth * s->nr * sizeof(double) );
  double * fpolar2_imag = malloc( s->nth * s->nr * sizeof(double) );


  for(j=0;j<s->nth;j++){
    imag[j] = 0.0;
    imag_b[j] = 0.0;
  }

  for(i=0;i<nr;i++){ 
    for(j=0;j<nth;j++){
      fpolar_real[i*nth+j] = 0.0;
      fpolar_imag[i*nth+j] = 0.0;

      fpolar2_real[i*nth+j] = 0.0;
      fpolar2_imag[i*nth+j] = 0.0;
    }
  }


  strcpy( fname, "wisdom.txt" );
  initialize_fftvars( &v, nth, fname );
  fftw_setup( &v );


  for(i=0;i<nr;i++){ 
    for(j=0;j<nth;j++){
      in1d[j] = polar[ i*nth+j];
      imag[j] = 0.0;
      in1d_b[j] = polar2[ i*nth+j];
      imag_b[j] = 0.0;
    }
  
    fftw( &v, in1d, imag, out1d, imag, nth, 1 );
    fftw( &v, in1d_b, imag_b, out1d_b, imag_b, nth, 1 );

    for(j=0;j<nth;j++){
      fpolar_real[ i*nth+j] = out1d[j];
      fpolar_imag[ i*nth+j] = imag[j];

      fpolar2_real[ i*nth+j] = out1d_b[j];
      fpolar2_imag[ i*nth+j] = imag_b[j];
    }

  }

  
  
  for(i=0;i<nr;i++){
    for(k=0;k<nr;k++){
    
      for(j=0;j<nth;j++){
	out1d[j] = fpolar_real[ i*nth+j];
	out1d_b[j] = fpolar2_real[ k*nth+j];
	imag[j] = fpolar_imag[ i*nth+j];
	imag_b[j] = fpolar2_imag[ k*nth+j];
      }
   
      for(j=0;j<nth;j++){
	tempr[j] = out1d[j] * out1d_b[j] + imag[j]*imag_b[j];
	tempi[j] = -out1d[j] * imag_b[j] + out1d_b[j]*imag[j];
      }

      fftw( &v, tempr, tempi, in1d, imag, nth, -1 );
      
      for(j=0;j<nth;j++){
	corr[ i*nth*nr+k*nth+j] = in1d[j]; //sqrt( in1d[j]*in1d[j] + imag[j]*imag[j]);
      }
      
    }
  }


  free(fpolar_real);
  free(fpolar_imag);
  free(fpolar2_real);
  free(fpolar2_imag);
  free(in1d);
  free(out1d);
  free(imag);
  free(in1d_b);
  free(out1d_b);
  free(imag_b);
  free(tempr);
  free(tempi);

}

