/*
 * qgrids.h      : 
 *      q sampling of the detector
 *
 * Andrew V. Martin March 2014
 *
 */


#ifndef QGRIDS_H
#define QGRIDS_H



#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "settings.h"


typedef struct{

  int nx;
  int nrho;
  double *qx, *qy, *qz;
  double *q2;
  double *qmod1d;
  
} qgrids;

/*
 *  allocate arrays for q grids
 *  calculate q-grids, q^2 for each detector pixel
 */
void calculate_qgrids( settings * s, qgrids * qg );

/*
 *  free allocated q grids
 */
void free_qgrids(qgrids * qg );


#endif
