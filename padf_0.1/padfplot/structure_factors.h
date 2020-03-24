/*
 * structure_factors.h      : 
 *      sample the structurefactors on the detector array
 *
 * Andrew V. Martin Jan 2013
 *
 */


#ifndef STUCT_FACT_H
#define STRUCT_FACT_H

#ifndef E_RAD
 #define CLASSICAL_E_RADIUS 2.8179403267*1e-15  // m
#endif


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf.h>
#include "bwxray.h"
#include "pdb.h"
#include "settings.h"
#include "io.h"
#include "geom_tools.h"

typedef struct{
  
  char code[3];
  int  Z;
  double * sf_array;
  double * sf_1d;

} element_label; 

typedef struct{

  int nel;
  element_label e[92];

} element_look_up;

typedef struct{

  int nel;
  element_label element[92];
  
} sf_list;

typedef struct{

  int nx;
  int nrho;
  double *qx, *qy, *qz;
  double *q2;
  double *qmod1d;
  
} qgrids;
 

/*
 *  Calculate a moment rho_lm(q)
 */
void calculate_moment( settings*s, pdbdata * pd, sf_list* sfl, 
		       qgrids * qg, int l, int m);


/*
 *  Calculate the diffraction pattern of the molecule
 */
void calculate_diffraction( settings* s, pdbdata * pd, sf_list* sfl, 
			    qgrids * qg);

/*
 *  Allocate the sf_list array
 *  Use the qgrids to calculate the scattering factors 
 *  for each element
 */
void calculate_structure_factors( settings * s, sf_data * Fdata, 
				  sf_list * sfl, qgrids * qg );

/*
 *  Free sf_list arrays
 */
void free_sf_list( sf_list * sfl);


/*
 *  allocate arrays for q grids
 *  calculate q-grids, q^2 for each detector pixel
 */
void calculate_qgrids( settings * s, qgrids * qg );

/*
 *  free allocated q grids
 */
void free_qgrids(qgrids * qg );


/*
 *  scan pdb data for unique elements
 *  store them in the structure factor list (sf_list)
 */
void pdb_unique_elements( pdbdata * pd , sf_list * sfl );

/*
 *  A list that assiates letter and numbers for elements
 */
void fill_element_look_up( element_look_up * el);

/*
 *  Evaluate a spherical harmonic at angles th and ph
 */

void sph_harmonic( int l, int m, double th, double ph, double* realout, double* imagout);

#endif
