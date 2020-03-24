/*
 * settings.c : part of padf
 *
 * A.V. Martin March 2015
 */


# include "settings.h"

void initialize_settings( settings * s)
{

  strcpy( s->outpath, "None");
  strcpy( s->input_dp_name, "None");
  strcpy( s->input_dp_name2, "None");
  strcpy( s->tag, "Tag");

  strcpy( s->logname, s->tag);
  strcat( s->logname, "padfcorr_log.txt");

  /*
   *  Diffaction parameters
   */
  s->cx = -1;  // negative value will default to centre of diffraction pattern array
  s->cy = -1;  // negative value will default to centre of diffraction pattern array
  s->nth = -1; // if negative set to ~ 2pi * s->nq
  s->nq = -1;  // if negative, set to half pixel side length of input array
  s->nr = -1;  // if negative, set to s->nq
  s->qmax = -1;   // If negative, it will be set to maximum detector q
  s->rmax = -1;   // If negative, it will be set to inverse detector pixel q width

  s->wavelength = 1e-10;
  s->pixel_width = 3.38e-5;
  s->detector_z = 1e-2;

  s->output_xsections = 0;
  s->Xflag = 0;

  s->nx = -1;
  s->ny = -1;
  s->pwx = -1;
  s->pwy = -1;
  

}


/*
 * Write settings to file
 */
void writeConfig(settings * s){

  /*
  printf("outpath = %s\n", s->outpath);
  printf("input diffraction pattern = %s\n", s->input_dp_name);
  printf("tag = %s\n", s->tag);
  printf("wavelength (m)                  = %g\n"  , s->wavelength  );
  printf("detector distance (m)           = %g\n"  , s->detector_z  );
  printf("pixel width (m)                 = %g\n"  , s->pixel_width  );
  printf("Pattern centre x                = %d\n"  , s->cx  );
  printf("Pattern centre y                = %d\n"  , s->cy  );
  printf("Samples q polar angle (theta)   = %d\n", s->nth );
  printf("Samples q radial                = %d\n", s->nq );
  printf("Samples r radial                = %d\n", s->nr );
  printf("qmax                            = %g\n", s->qmax );
  printf("rmax                            = %g\n", s->rmax );
  */

  FILE * outfile = fopen( s->logname, "a" );

  if (outfile != NULL) {  
    fprintf(outfile,"outpath = %s\n", s->outpath);
    fprintf(outfile,"input diffraction pattern = %s\n", s->input_dp_name);
    fprintf(outfile,"input diffraction pattern 2 = %s\n", s->input_dp_name2);
    fprintf(outfile,"tag = %s\n", s->tag);
    fprintf(outfile,"wavelength (m)                  = %g\n"  , s->wavelength  );
    fprintf(outfile,"detector distance (m)           = %g\n"  , s->detector_z  );
    fprintf(outfile,"pixel width (m)                 = %g\n"  , s->pixel_width  );
    fprintf(outfile,"Pattern centre x                = %d\n"  , s->cx  );
    fprintf(outfile,"Pattern centre y                = %d\n"  , s->cy  );
    fprintf(outfile,"Samples q polar angle (theta)   = %d\n", s->nth );
    fprintf(outfile,"Samples q radial                = %d\n", s->nq );
    fprintf(outfile,"Samples r radial                = %d\n", s->nr );
    fprintf(outfile,"qmax                            = %g\n", s->qmax );
    fprintf(outfile,"rmax                            = %g\n", s->rmax );
    fprintf(outfile,"output_xsections                = %d\n", s->output_xsections );
    fprintf(outfile,"cross-correlation flag (Xflag)  = %d\n", s->Xflag );
    fprintf(outfile,"No. pixels x                    = %d\n", s->nx );
    fprintf(outfile,"No. pixels y                    = %d\n", s->ny );
    fprintf(outfile,"pixel width x (pwx)             = %g\n", s->pwx );
    fprintf(outfile,"pixel width y (pwy)             = %g\n", s->pwy );
  } 
  
  fclose( outfile );
}


/*
 *  Read a config file
 */

void readConfig(char * fname, settings * s){
  
  FILE * file = fopen ( fname, "r" );
  char line[1024];
  int i=0;
  line[0] = '\0';

  if (file != NULL){
    while( fgets(line, sizeof line, file) != NULL) {

      if(line[0] != '#'){
	parseLine(line, s);
      }
    }
   
    fclose(file);
  } 
  else {
    printf("Config file was not found.\n");
  }


  /*
   *  parameters and arrays for moments
   */
  strcpy( s->logname, s->outpath );
  strcat( s->logname, s->tag );
  strcat( s->logname, "_log.txt");
  file = fopen( s->logname,   "w");
  if (file == NULL){
    printf("log file could not be opened - check the directory given in the config file\n");
    printf("%s\n", s->logname);
  } else {
    fclose( file );
  }


}


void parseLine(char * line, settings * s){

  
  char label[1024];
  char value[1024];

  char a = 'a';
  char *  pa = &a;
  char ** endptr = &pa;

  FILE * fp;


  label[0] = '\0';
  value[0] = '\0';

  breakLine(line, label, value);
  

  if (strcmp(label,"wavelength") == 0) {
    s->wavelength = atof(value); 
  } 
  else if (strcmp(label,"outpath") == 0) {
    strcpy(s->outpath, value);
    if ( (s->outpath[strlen(s->outpath)-1] != '/')
	 || (s->outpath[strlen(s->outpath)-1] != '\\') ){
      strcat( s->outpath, "/");
    }
  }
  else if (strcmp(label,"tag") == 0) {
    strcpy(s->tag, value);
  }
  else if (strcmp(label,"input") == 0) {
    strcpy(s->input_dp_name, value);
    pa = value;
    fp = fopen( pa,   "r");
    if (fp == NULL){
      printf("diffraction pattern file could not be opened - check the directory/name given in the config file\n");
      printf("%s\n", s->input_dp_name);
    } else {
      fclose( fp );
    }
  }
  else if (strcmp(label,"inputb") == 0) {
    strcpy(s->input_dp_name2, value);
    pa = value;
    fp = fopen( pa,   "r");
    if (fp == NULL){
      printf("diffraction pattern file 2 could not be opened - check the directory/name given in the config file\n");
      printf("%s\n", s->input_dp_name2);
    } else {
      fclose( fp );
    }
  }
  else if (strcmp(label,"detector_z") == 0) {
    s->detector_z = atof(value); 
  }
  else if (strcmp(label,"pixel_width") == 0) {
    s->pixel_width = atof(value); 
  }
  else if (strcmp(label,"cx") == 0) {
    s->cx = atoi(value); 
  }
  else if (strcmp(label,"cy") == 0) {
    s->cy = atoi(value); 
  }
  else if (strcmp(label,"nth") == 0) {
    s->nth = atoi(value); 
  }
  else if (strcmp(label,"nq") == 0) {
    s->nq = atoi(value); 
  }
  else if (strcmp(label,"nr") == 0) {
    s->nr = atoi(value); 
  }
  else if (strcmp(label,"qmax") == 0) {
    s->qmax = atof(value); 
  }
  else if (strcmp(label,"rmax") == 0) {
    s->rmax = atof(value); 
  }
  else if (strcmp(label,"output_xsections") == 0) {
    s->output_xsections = atoi(value); 
  }
  else if (strcmp(label,"xflag") == 0) {
    s->Xflag = atoi(value); 
  }
  else if (strcmp(label,"nx") == 0) {
    s->nx = atoi(value); 
  }
  else if (strcmp(label,"ny") == 0) {
    s->ny = atoi(value); 
  }
  else if (strcmp(label,"pwx") == 0) {
    s->pwx = atof(value); 
  }
  else if (strcmp(label,"pwy") == 0) {
    s->pwy = atof(value); 
  }
  else if (strcmp(label,"none") == 0) {
    ;
  }
  else
  {
    printf("input line not parsed : %s \n", line);
  }

}


void breakLine(char * line, char * label, char * value){
  char * token;
  token = "\0";

  token = strtok(line," =:\n");
  if ( token != NULL ){
    convertToLower(token);
    strcpy( label, token );
  } else {
    strcpy( label, "none" );
  }

  if( token != NULL ){
    token = strtok(NULL," =:\n");
    strcpy( value, token );
  } else {
    strcpy( value, "none" );
  }


  //printf("Label %s  Value %s",label, value);
  
}

void convertToLower(char * str)
{
  int len = 0;
  int i=0;
  len = strlen("erer");
 
  for (i=0; i<len; i++)
    {
      str[i] = tolower(str[i]);
    }

}


void parse_config_name( int argc, char * argv[], char * cname ){

  if ( argc >= 2 ){

    strcpy( cname, argv[1] );

  }

}


char *trimwhitespace(char *str)
{
  char *end;

  // Trim leading space
  while(isspace(*str)) str++;

  if(*str == 0)  // All spaces?
    return str;

  // Trim trailing space
  end = str + strlen(str) - 1;
  while(end > str && isspace(*end)) end--;

  // Write new null terminator
  *(end+1) = 0;

  return str;
}


void default_settings( settings * s, int n )
{

  /*
   *  Diffraction parameters
   */
  if (s->nx==-1) {
    s->nx = sqrt(n);
  }

  if (s->ny==-1) {
    s->ny = sqrt(n);
  }

  int nx = s->nx;
  int ny = s->ny;

  if( s->cx < 0){
    s->cx = nx/2;
  }

  if( s->cy < 0){
    s->cy = ny/2;
  }


  if(s->nq < 0){
    s->nq = nx/2;
  }

  if(s->nr < 0){
    s->nr = s->nq;
  }

  if (s->nth < 0){
    s->nth = (int) 2 * 3.14159 * (s->nr/2);
  }
 
  if(s->qmax < 0){
    double thmax = atan ( (nx/2) * s->pixel_width / s->detector_z );
    s->qmax = (2/s->wavelength) * sin( thmax/2.0 );
  }

  if(s->rmax < 0){
    s->rmax = 0.25*s->nr/s->qmax;
  }

  if(s->pwx < 0.0){
    s->pwx = s->pixel_width;
  }

  if(s->pwy < 0.0){
    s->pwy = s->pixel_width;
  }

}
