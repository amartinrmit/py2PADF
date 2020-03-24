
#include "qgrids.h"

/*
 *  calculate q-grids, q^2 for each detector pixel
 */
void calculate_qgrids( settings * s, qgrids * qg )
{
  
  int i, j;
  double *x, *y, *r;
  double length;

  /*
   *  set qgrid array size
   */
  qg->nx   = s->nx;
  qg->nrho = s->nx/2; //s->nrho;

  
  /*
   * allocate arrays
   */
  qg->qx = malloc( qg->nx*qg->nx *sizeof(double) );
  qg->qy = malloc( qg->nx*qg->nx *sizeof(double) );
  qg->qz = malloc( qg->nx*qg->nx *sizeof(double) );
  qg->q2 = malloc( qg->nx*qg->nx *sizeof(double) );
  qg->qmod1d = malloc( qg->nrho  *sizeof(double) );

  x = malloc( qg->nx*qg->nx *sizeof(double) );
  y = malloc( qg->nx*qg->nx *sizeof(double) );
  r = malloc( qg->nrho      *sizeof(double) );


  /*
   * get coordinates for each for each detector pixel
   */
  for (i=0;i<qg->nx;i++){
    for (j=0;j<qg->nx;j++){
      x[i*qg->nx+j] = (j * s->pixel_width) - s->cx;
      y[i*qg->nx+j] = (i * s->pixel_width) - s->cy;
    }
  }

  for (j=0;j<qg->nrho;j++){
    r[j] = (j * s->pixel_width);
  }


  /*
   * calculate the q grids
   */
  for (i=0;i<qg->nx;i++){
    for (j=0;j<qg->nx;j++){
      length = sqrt( x[i*qg->nx+j]*x[i*qg->nx+j]
	+ y[i*qg->nx+j]*y[i*qg->nx+j]
	+ s->detector_z*s->detector_z );

      qg->qx[i*qg->nx+j] = x[i*qg->nx+j] / ( length*s->wavelength*1e10 );
      qg->qy[i*qg->nx+j] = y[i*qg->nx+j] / ( length*s->wavelength*1e10 );    
      qg->qz[i*qg->nx+j] = - (1./(s->wavelength*1e10))
	+ s->detector_z / ( length*s->wavelength*1e10 );
      
      qg->q2[i*qg->nx+j] = sqrt( qg->qx[i*qg->nx+j]*qg->qx[i*qg->nx+j]
	+ qg->qy[i*qg->nx+j]*qg->qy[i*qg->nx+j]
	+ qg->qz[i*qg->nx+j]*qg->qz[i*qg->nx+j] );

    }
  }

  for (j=0;j<qg->nrho;j++){
      length = sqrt( r[j]*r[j] + s->detector_z*s->detector_z );

      qg->qmod1d[j] = r[j] / ( length*s->wavelength*1e10 );
  }

 
  /*
   *  Output qgrids
   */
  /*
  char outname[1024];
  strcpy( outname, s->outpath );
  strcat( outname, "/q2.dbin");
  write_dbin( outname, qg->q2, s->nx*s->nx );
  */


  /*
   * free arrays
   */
  free(x);
  free(y);

}


/*
 *  free allocated q grids
 */
void free_qgrids(qgrids * qg )
{
  free(qg->qx);
  free(qg->qy);
  free(qg->qz);
  free(qg->q2);
  free(qg->qmod1d);
}
