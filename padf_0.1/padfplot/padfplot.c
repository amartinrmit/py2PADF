
/*
 *  padfplot.c 
 *  A.V. Martin  March 2015
 *
 *  Calculate the pair-angle distribution function from a real-space density matrices
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
#include "plotting.h"

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
   * allocate and read in density matrices
   */
  blarr_t blar;
  allocate_blarr_t( &blar, s.nl, s.nr );
  // read in the blarr functions...
  for(i=0;i<s.nl;i+=2){
    set_fname_numbered( &iov, "blrr", ".dbin", i, outname );
    read_dbin( outname, blar.l[i].data, s.nr*s.nr );
  }
  printf("DEBUG s.nr %d\n", s.nr);

  double * thplane = malloc( s.nr*s.nr*sizeof(double) );
  double * rplane = malloc( s.nr*s.nthq*sizeof(double) );
  if( s.section == 1 ){

    for(i=0;i<s.nr*s.nr;i++){ thplane[i] = 0.0; }

    // calculate a plane of constant theta
    calculate_thplane( &blar, s.nthq, s.theta, s.rmax, s.nl, thplane);

    // output the plane
    strcpy( outname, s.outpath );
    strcat( outname, "/" );
    strcat( outname, s.tag );
    strcat( outname, "_thplane.dbin" );
    write_dbin( outname, thplane, s.nr*s.nr );
   

  } else if (s.section == 2 ){


    for(i=0;i<s.nr*s.nthq;i++){ rplane[i] = 0.0; }
    printf("DEBUG rplane r %g\n", s.r);
    printf("DEBUG rmax r %g\n", s.rmax);
    printf("DEBUG rplane nl %d\n", s.nl);
    printf("DEBUG rplane nthq %d\n", s.nthq);

    // calculate a plane of constant r
    calculate_rplane( &blar, s.nthq, s.r, s.rmax, s.nl, rplane);

    // output the plane
    strcpy( outname, s.outpath );
    strcat( outname, "/" );
    strcat( outname, s.tag );
    strcat( outname, "_rplane.dbin" );
    write_dbin( outname, rplane, s.nr*s.nthq );
    // free(rplane);


  } else if (s.section == 3 ){

    //    double * rplane = malloc( s.nr*s.nthq*sizeof(double) );

    // calculate a plane of constant r
    calculate_r_eq_r_plane( &blar, s.nthq, s.rmax, s.nl, rplane);

    // output the plane
    strcpy( outname, s.outpath );
    strcat( outname, "/" );
    strcat( outname, s.tag );
    strcat( outname, "_rrplane.dbin" );
    write_dbin( outname, rplane, s.nr*s.nthq );
    // free(rplane);



  } else if (s.section == 4 ){

    double * thline = malloc( s.nthq*sizeof(double) );

    // calculate a line out as a function of theta
    calculate_thline( &blar, s.nthq, s.r, s.r2, s.rmax, s.nl, thline);

    // output the plane
    strcpy( outname, s.outpath );
    strcat( outname, "/" );
    strcat( outname, s.tag );
    strcat( outname, "_thline.dbin" );
    write_dbin( outname, thline, s.nthq );
    free(thline);

    
  } else if (s.section == 5 ){

    double * rline = malloc( s.nr*sizeof(double) );

    // calculate a line out as a function of theta
    calculate_rline( &blar, s.theta, s.nthq, s.r, s.rmax, s.nl, rline);
    
    // output the plane
    strcpy( outname, s.outpath );
    strcat( outname, "/" );
    strcat( outname, s.tag );
    strcat( outname, "_rline.dbin" );
    write_dbin( outname, rline, s.nr );
    free(rline);

  }


  /*
   * free memory
   */
  clean_up_blarr_t( &blar );
  free(thplane);
  free(rplane);
}
