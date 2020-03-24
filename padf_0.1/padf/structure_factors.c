/*
 * structure_factors.c      : 
 *      sample the structurefactors on the detector array
 *
 * Andrew V. Martin Jan 2013
 *
 */


#include "structure_factors.h"


/*
 *  Calculate a moment rho_lm(q)
 */
void calculate_moment( settings*s, pdbdata * pd, sf_list* sfl, 
		       qgrids * qg, int l, int m)
{
  int i=0;

  //sph_harmonic( int l, int m, double th, double ph, double* realout, double* imagout);
  // gsl_sf_bessel_jl (int l, double x);

  int ia, ip, ie;
  int Zindex;
  double dotp;
  double vect[3];
  double x, y, cosTheta, sa;
  double R, th, ph;
  double sphH_real, sphH_imag, sphB;


  /*
   *  initialize the output arrays
   */
  for (ip=0;ip<s->nrho;ip++){  
    s->rho_real[ip] = 0.0;
    s->rho_imag[ip] = 0.0;
  }    

  double dq = (1.0*s->wavelength)*s->pixel_width/s->detector_z;

  /*
   * Sum over all the atoms
   */
  for (ia=0;ia<pd->natoms;ia++){
    
    /*
     *  Progress update to screen
     */
    if (ia % 100 == 0){
      printf( "Calculating moment... Atoms added %d / %d\n", ia, pd->natoms);
    }


    /*
     *  Identify the atom species
     */
    for (ie=0;ie<sfl->nel;ie++){
    	if (strcmp( pd->adata[ia].element,  
		    sfl->element[ie].code ) == 0){
	  Zindex = ie;
	  break;
	}     
    }


    /*
     *  Loop over all pixels
     */
    for (ip=0;ip<s->nrho;ip++){
     
      vect[0] = pd->adata[ia].x;
      vect[1] = pd->adata[ia].y;
      vect[2] = pd->adata[ia].z;
     
      /*
       *  If specified, rotate the molecule
       */
      if (s->angle != 0.0){
	rotate_vector( vect, s->axis, s->angle, vect );
      }
      
      /*
       *  Atom location in spherical coordinates
       */
      R = sqrt(vect[0]*vect[0] + vect[1]*vect[1] + vect[2]*vect[2]);

      if (abs(R) > 0.0 ){
	th = acos( vect[2] / R );
      } else {
	th = 0.0;
      }

      if ( vect[0] > 0.0){
	ph = atan( vect[1] / vect[0] ) + PI/2;
      } else if ( vect[0] < 0.0){
	ph = atan( vect[1] / vect[0] ) - PI/2;
      } else {
	ph = 0.0;
      }
   
      /*
       * Calculate the moment
       */
      sph_harmonic( l, m, th, ph, &sphH_real, &sphH_imag);
      sphB =  gsl_sf_bessel_jl (l, qg->qmod1d[ip]*R);

      s->rho_real[ip] += sfl->element[Zindex].sf_1d[ip] * sphB * sphH_real  ;
      s->rho_imag[ip] += sfl->element[Zindex].sf_1d[ip] * sphB * sphH_imag  ;

    } //endfor ip   

  } //endfor ia
  



}    



/*
 *  Calculate the diffraction pattern of the molecule
 */
void calculate_diffraction( settings* s, pdbdata * pd, sf_list* sfl, 
			    qgrids * qg)
{

  int ia, ip, ie;
  int Zindex;
  double dotp;
  double vect[3];
  double x, y, cosTheta, sa;

  /*
   *  initialize the output arrays
   */
  for (ip=0;ip<s->nx*s->nx;ip++){  
    s->sf_real[ip] = 0.0;
    s->sf_imag[ip] = 0.0;
  }    



  /*
   * Sum over all the atoms
   */
  for (ia=0;ia<pd->natoms;ia++){
    
    /*
     *  Progress update to screen
     */
    if (ia % 100 == 0){
      printf( "Atoms added %d / %d\n", ia, pd->natoms);
    }


    /*
     *  Identify the atom species
     */
    for (ie=0;ie<sfl->nel;ie++){
    	if (strcmp( pd->adata[ia].element,  
		    sfl->element[ie].code ) == 0){
	  Zindex = ie;
	  break;
	}     
    }


    /*
     *  Loop over all pixels
     */
    for (ip=0;ip<s->nx*s->nx;ip++){
     
      vect[0] = pd->adata[ia].x;
      vect[1] = pd->adata[ia].y;
      vect[2] = pd->adata[ia].z;
     
      /*
       *  If specified, rotate the molecule
       */
      if (s->angle != 0.0){
	rotate_vector( vect, s->axis, s->angle, vect );
      }

      /*
       *  Calculate the dot product with pixel direction
       */
      dotp = vect[0] * qg->qx[ip]
	+ vect[1] * qg->qy[ip]
	+ vect[2] * qg->qz[ip];
	
     
      /*
       * Calculate the structure factor, real and imaginary parts
       */
      s->sf_real[ip] += 
	sfl->element[Zindex].sf_array[ip]*cos( 2.0 * PI * dotp  );
      s->sf_imag[ip] += 
	sfl->element[Zindex].sf_array[ip]*sin( 2.0 * PI * dotp  );
    }    

  } //endfor ia
  
  /*
   *  Calculate the diffraction pattern
   */
  for (ip=0;ip<s->nx*s->nx;ip++){  

     /*
      *  Calculate the solid angle effect
      */
    x = ((ip%s->nx) * s->pixel_width) - s->centre_x;
    y = ((ip/s->nx) * s->pixel_width) - s->centre_y;
    cosTheta = cos( atan(sqrt(x*x + y*y) / s->detector_z) );
    sa = pow( s->pixel_width/s->detector_z, 2 ) * cosTheta;
    //DEBUG**************
    //printf("ip x y sa %d %g %g %g %g\n", ip, x, y, cosTheta, sa );

    /*
     * Evaluate intensity
     */
    s->intensity[ip] = sa*(s->sf_real[ip]*s->sf_real[ip]
			+ s->sf_imag[ip]*s->sf_imag[ip]);
  }    

  

}



/*
 *  Allocate the sf_list array
 *  Use the qgrids to calculate the scattering factors 
 *  for each element
 */
void calculate_structure_factors( settings * s, sf_data * Fdata, 
				  sf_list * sfl, qgrids * qg )
{

  int ix, iy, ie;
  int Z, ic;
  double c1, c2;


  /*
   *  Loop over all the unique elements
   */
  for (ie=0; ie<sfl->nel; ie++){
    
    Z = sfl->element[ie].Z-1;

    /*
     *  Allocate the sfl array
     */
    sfl->element[ie].sf_array = malloc( qg->nx*qg->nx*sizeof(double) );
    sfl->element[ie].sf_1d    = malloc( qg->nrho     *sizeof(double) );


    /*
     *  Loop over all pixels on the detector
     */
    for (ix=0; ix<qg->nx; ix++) {
      for (iy=0; iy<qg->nx; iy++) {
	
	/*
	 *  calculate the scattering factor for this pixel
	 */
	sfl->element[ie].sf_array[ix*qg->nx+iy] = 0.0;

	for (ic=0;ic<5;ic++){
	  c1 = Fdata->element[Z].coeff[ic];
	  c2 = Fdata->element[Z].coeff[ic+5];
	  sfl->element[ie].sf_array[ix*qg->nx+iy] += 
	    c1*exp(-c2*qg->q2[ix*qg->nx+iy]/4.0 );
	}
	  sfl->element[ie].sf_array[ix*qg->nx+iy] += 
	    Fdata->element[Z].coeff[10];

      } //end for iy
    }  // end for ix




    /*
     *  Loop over all pixels in 1D radial line
     */
    for (iy=0; iy<qg->nrho; iy++) {
    
      /*
       *  calculate the scattering factor for this pixel
       */
      sfl->element[ie].sf_1d[iy] = 0.0;
    
      for (ic=0;ic<5;ic++){
	c1 = Fdata->element[Z].coeff[ic];
	c2 = Fdata->element[Z].coeff[ic+5];
	sfl->element[ie].sf_1d[iy] += 
	  c1*exp(-c2*pow(qg->qmod1d[iy],2)/4.0 );
      }
      sfl->element[ie].sf_1d[iy] += 
	Fdata->element[Z].coeff[10];

    } //end for iy

  } //end for ie

  /*
   *  Output qgrids
   */
  /*
  char outname[1024];
  strcpy( outname, s->outpath );
  strcat( outname, "/scattering_factor.dbin");
  write_dbin( outname, sfl->element[1].sf_array, s->nx*s->nx );
  */



}


/*
 *  Free sf_list arrays
 */
void free_sf_list( sf_list * sfl)
{

  int ie;

  for (ie=0; ie<sfl->nel; ie++){
    free( sfl->element[ie].sf_array);
    free( sfl->element[ie].sf_1d);
  }

}


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
  qg->nrho = s->nrho;

  
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
      x[i*qg->nx+j] = (j * s->pixel_width) - s->centre_x;
      y[i*qg->nx+j] = (i * s->pixel_width) - s->centre_y;
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




/*
 *  scan pdb data for unique elements
 */
void pdb_unique_elements( pdbdata * pd ,  sf_list * sfl )
{
  int i,j;
  int unique = 1;

  sfl->nel = 0;

  element_look_up elist;
  fill_element_look_up( &elist );

   
  for (i=0;i<pd->natoms;i++) {

    unique = 1;
    /*
     *  Search for the element 'i' in the current list
     */
    for (j=0;j<sfl->nel;j++) {
      if (strcmp( pd->adata[i].element, sfl->element[j].code ) == 0){
	unique = 0;
      }
    }

    /*
     *  If a new element has been identified, add it to the list
     */
    if (unique == 1){
      strcpy( sfl->element[sfl->nel].code, pd->adata[i].element );

      /*
       *  Find the numeric label in the elist
       */
      int found = 0;
      for (j=0;j<92;j++){
	
	if (strcmp( elist.e[j].code,  
		    sfl->element[sfl->nel].code ) == 0){
	  sfl->element[sfl->nel].Z = elist.e[j].Z;
	  found = 1;
	  break;
	}
      }
      if (found == 0){
	printf("ERROR: pdb code %s could be found in the element list\n", 
	       sfl->element[sfl->nel].code);
	sfl->element[sfl->nel].Z = 0;
      }

      /*
       * a new element was found, so increment list size
       */
      sfl->nel ++;
      
    } //endif unique

  } // end loop: pd.natoms
 
  printf("Number of unique elements = %d\n", sfl->nel);

}





/*
 *  A list that assiates letter and numbers for elements
 */
void fill_element_look_up( element_look_up * el)
{
  
  int i;
  char names[] = " H, He,"
    "Li,Be, B, C, N, O, F,Ne,"
    "Na,Mg,Al,Si, P, S,Cl,Ar,"
    " K,Ca,Sc,Ti, V,Cr,Mn,Fe,Co,Ni,Cu,Zn,Ga,Ge,As,Se,Br,Kr,"
    "Rb,Sr, Y,Zr,Nb,Mo,Tc,Ru,Rh,Pd,Ag,Cd,In,Sn,Sb,Te, I,Xe,"
    "Cs,Ba,"
    "La,Ce,Pr,Nd,Pm,Sm,Eu,Gd,Tb,Dy,Ho,Er,Tm,Yb,"
    "Lu,Hf,Ta, W,Re,Os,Ir,Pt,Au,Hg,Tl,Pb,Bi,Po,At,Rn,"
    "Fr,Ra,"
    "Ac,Th,Pa, U";

  char * token;
  token = "\0";
  
  strtok( names, "," );

  for (i=0;i<92;i++){
    strcpy( el->e[i].code, token);
    el->e[i].Z = i+1;
    el->e[i].sf_array = NULL;
    token = strtok(NULL,",\n"); 
  }

  el->nel = 92;

}

/*
 *  Evaluate a spherical harmonic at angles th and ph
 */

void sph_harmonic( int l, int m, double th, double ph, double* realout, double* imagout)
{

  double assocL; 
  gsl_complex output, temp;

  int mabs = labs(m);
 
  assocL = gsl_sf_legendre_sphPlm(l, mabs, cos(th) );
 
  if (m < 0){
    assocL *= pow(-1.0, mabs ) ;
  }
  
  *realout = assocL*cos(m*ph );  
  *imagout = assocL*sin(m*ph );  
 
}
