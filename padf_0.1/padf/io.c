

#include "io.h"


void write_dbin( char * fname, double * array, int n){

  FILE * fp;

  fp=fopen(fname, "wb");
  
  if (fp != NULL){
    fwrite(array, sizeof(double), n, fp);
    fclose(fp);
  }

}

void write_int_to_dbin( char * fname, 
			int * array, int n){

  FILE * fp;
  double * darray;
  darray = malloc( n*sizeof(double));

  int i;
  for (i=0;i<n;i++){
    darray[i] = 1.0*array[i];
  }


  fp=fopen(fname, "wb");
  
  if (fp != NULL){
    fwrite( darray, sizeof(double), n, fp);
    fclose(fp);
  }

  free(darray);
}


void read_dbin( char * fname, double * array, int n){

  FILE * fp;

  fp=fopen(fname, "rb");
  
  fread(array, sizeof(double), n, fp);

  fclose(fp);
}


long number_of_pixels_in_file( char * fname ){
  
  FILE * f;
  long fsize;
  f = fopen( fname, "rb");
  fseek( f, 0, SEEK_END);
  fsize = ftell(f) / sizeof(double) ;
  fclose(f);

  return fsize;
}

void allocate_and_read_binary_image( char * fname, long * n, double ** image){

  FILE * f;
  if(*image != NULL){ 
    free(*image);
  }

  *n = number_of_pixels_in_file( fname );
  *image = malloc( *n * sizeof(double) );
  
  f = fopen( fname, "rb");
  fread( *image, sizeof(double), *n, f );
  fclose(f);
}

void read_binary_image( char * fname, long n, double * image){

  FILE * f;
  
  f = fopen( fname, "rb");
  fread( image, sizeof(double), n, f );
  fclose(f);
}


void set_io_vars_t( io_vars_t * iovars, char * outpath, char * tag ){

  strcpy( iovars->outpath, outpath );
  strcpy( iovars->tag, tag );

}

void set_fname( io_vars_t * iovars, char * label, char * fext, char * fname ){

  strcpy( fname, iovars->outpath);
  strcat( fname, iovars->tag    );
  strcat( fname, label          );
  strcat( fname, fext           );

}


void set_fname_numbered( io_vars_t * iovars, char * label, char * fext, int i, char * fname  ){

  char str[10];
  sprintf(str, "_%d", i );

  strcpy( fname, iovars->outpath);
  strcat( fname, iovars->tag    );
  strcat( fname, label          );
  strcat( fname, str            );
  strcat( fname, fext           );

}
