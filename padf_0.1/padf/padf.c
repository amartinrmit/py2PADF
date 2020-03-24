
/*
 *  padf.c 
 *  A.V. Martin  March 2015
 *
 *  Calculate the pair-angle distribution function from a diffraction pattern
 *
 */


#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include "settings.h"
#include "io.h"
#include "bll.h"

int main(int argc, char * argv[]){

  int i;
  settings s;
  char outname[1024];
  char configname[1024];
  strcpy( configname, "config.txt");
  parse_config_name( argc, argv, configname );

  //srand(time(NULL));

  /*
   *  Initialize file names, experiment parameters
   */
  initialize_settings( &s );
  readConfig( configname, &s);

  io_vars_t iov;
  iov.outpath[0] = '\0';
  //set_io_vars_t( &iov, s.outpath, s.tag );
  strcpy( iov.outpath, s.outpath);
  strcpy( iov.tag, s.tag);
  printf("DEBUG iov %s\n", iov.outpath);
  printf("DEBUG iov %s\n", iov.tag);
  printf("DEBUG settings %s %s\n", s.outpath, s.tag);

  /*
   * allocate and read in correlation function
   */
  correlation_t corr;
  allocate_correlation_t( &corr, s.nq, s.nthq );
  correlation_t corrsig;
  allocate_correlation_t( &corrsig, s.nq, s.nthq );
 
  read_dbin( s.correlationfile, corr.data, corr.npix );

  if (s.nfilterflag == 1 ){
    read_dbin( s.correlation_sigma_file, corrsig.data, corrsig.npix );
  } else {
    for (i=0;i<corrsig.npix;i++){ corrsig.data[i] = 0.0; }
  }

  /*
   * Allocate and read gfilter file
   */
  double * gfilter = malloc( s.nq * sizeof(double) );
  for (i=0;i<s.nq;i++){ gfilter[i] = 0.0; }
  if (s.gfilter_flag == 1){
    read_dbin( s.gfilterfile, gfilter, s.nq );
  }

  /*
   * allocate and read in noise data if required
   */
  noisedata_t nd;
  allocate_noisedata_t( &nd, s.nq );
  if (s.noise_estimation_flag==1){
    printf("nd.n %d\n", nd.n);
    read_dbin( s.qnoisefile, nd.sigma, nd.n );
    for (i=0;i<5;i++){ printf("noisedata %g\n", nd.sigma[i]); }
  } else {
    for (i=0;i<s.nq;i++){ nd.sigma[i] = 0.0; }
  }


  /*
   *  Print settings to screen and write to log
   */
  writeConfig( &s );

  /*
   *  Allocate padf variables
   */
  correlation_t padf;
  allocate_correlation_t( &padf, s.nr, s.nthq );
  correlation_t padf_noise;
  if (s.noise_estimation_flag == 1){
    allocate_correlation_t( &padf_noise, s.nr, s.nthq );
  } else {
    allocate_correlation_t( &padf_noise, 1, 1 );
    padf_noise.data[0] = -1;
  }
 
  /*
   * Calculate the padf
   */
  padf_calc( &corr, &corrsig, &padf, &padf_noise, s.rmax, s.qmax, s.wl, s.nl, s.nlmin, s.blflag, 
	     s.blfilter, s.nfilterflag, &iov, &nd, s.use_rl_filter, gfilter, s.gfilter_flag,
	     s.noise_estimation_flag);
   
  /*
   *  Output correlation
   */
  if(s.output_padf==1){
    strcpy( outname, s.outpath );
    strcat( outname, "/" );
    strcat( outname, s.tag );
    strcat( outname, "_padf.dbin" );
    write_dbin( outname, padf.data, padf.npix );
  }
  
  if (s.noise_estimation_flag==1){
    strcpy( outname, s.outpath );
    strcat( outname, "/" );
    strcat( outname, s.tag );
    strcat( outname, "_padf_noise.dbin" );
    write_dbin( outname, padf_noise.data, padf_noise.npix );
  }

  /*
   * free memory
   */
  clean_up_correlation_t( &padf );
  clean_up_correlation_t( &padf_noise );
  clean_up_correlation_t( &corr );
  clean_up_correlation_t( &corrsig );
  clean_up_noisedata_t( &nd );
  free(gfilter);

}
