
/*
 *  padfcorr.c 
 *  A.V. Martin  March 2015
 *
 *  Calculate the angular correlation function from a diffraction pattern
 *
 */


#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include "settings.h"
#include "io.h"
#include "qgrids.h"
#include "interpolation.h"
#include "correlation.h"

int main(int argc, char * argv[]){

  int i,j;
  settings s;
  char outname[1024];
  char configname[1024];
  strcpy( configname, "config.txt");
  parse_config_name( argc, argv, configname );
  //printf( "%d %s %s\n", argc, argv[1], configname );

  //srand(time(NULL));

  /*
   *  Initialize file names, experiment parameters
   */
  initialize_settings( &s );
  readConfig( configname, &s);

  /*
   * allocate and read in diffraction pattern image
   */
  double * dpin = NULL;
  long n;
  int nx;
  allocate_and_read_binary_image( s.input_dp_name, &n, &dpin );
  s.nx = sqrt(n);
  //  printf("DEBUG nx %d\n", s.nx);
  default_settings( &s, s.nx );

  double * dpin2 = NULL;
  allocate_and_read_binary_image( s.input_dp_name2, &n, &dpin2 );

  /*
   *  Print settings to screen and write to log
   */
  writeConfig( &s );

  /*
   *  Allocate polar and correlation variables
   */
  double * polar = malloc( s.nr * s.nth * sizeof(double) ) ;
  double * polar2 = malloc( s.nr * s.nth * sizeof(double) ) ;
  double * corr  = malloc( s.nr * s.nr * s.nth * sizeof(double) ) ;

  /*
   * Create polar plot
   */
  polar_plot( &s, dpin, polar );
  polar_plot( &s, dpin2, polar2 );

  /*
   *  Calculate correlation 
   */
  correlation_rtheta_faster( &s, polar, polar2, corr );
  
  /*
   *  Output correlation
   */
  strcpy( outname, s.outpath );
  strcat( outname, "/" );
  strcat( outname, s.tag );
  strcat( outname, "_correlation.dbin" );
  //printf("DEBUG outname %s\n", outname);
  //clock_t start = clock();
  write_dbin( outname, corr, s.nr*s.nr*s.nth );
  //clock_t end = clock();
  //printf("DEBUG writing correlation took %g seconds\n", (double) (end-start)/CLOCKS_PER_SEC );


  /*
   * Output cross-sections
   */
  if (s.output_xsections==1){
    //printf("DEBUG output xsections\n");
    double * qqsection = malloc( s.nr*s.nr*sizeof(double) );
    for (i=0;i<s.nr;i++){
      for (j=0;j<s.nr;j++){
	qqsection[i*s.nr+j] = corr[i*s.nr*s.nth+j*s.nth+0];
      }
    }
    strcpy( outname, s.outpath );
    strcat( outname, "/" );
    strcat( outname, s.tag );
    strcat( outname, "_correlation_theta0_section.dbin" );
    write_dbin( outname, qqsection, s.nr*s.nr );
    free(qqsection);

    double * qthsection = malloc( s.nr*s.nth*sizeof(double) );
    for (i=0;i<s.nr;i++){
      for (j=0;j<s.nth;j++){
	qthsection[i*s.nth+j] = corr[i*s.nr*s.nth+i*s.nth+j];
      }
    }
    strcpy( outname, s.outpath );
    strcat( outname, "/" );
    strcat( outname, s.tag );
    strcat( outname, "_correlation_q_eq_q_section.dbin" );
    write_dbin( outname, qthsection, s.nr*s.nth );
    free(qthsection);

  }



  /*
   * free memory
   */
  free( dpin );
  free( dpin2 );
  free( polar );
  free( polar2 );
  free( corr );

}
