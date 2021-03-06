/*
 * settings.c : part of padf
 *
 * A.V. Martin March 2015
 */


# include "settings.h"

void initialize_settings( settings * s)
{

  strcpy( s->outpath, "None");
  strcpy( s->correlationfile, "None");
  strcpy( s->correlation_sigma_file, "None");
  strcpy( s->qfile, "None");
  strcpy( s->tag, "Tag");

  strcpy( s->logname, s->tag);
  strcat( s->logname, "_log.txt");

  /*
   *  Diffraction parameters
   */
  s->nthq = -1; 
  s->nq = -1;  
  s->nr = -1;  
  s->qmax = -1;
  s->rmax = -1;   
  s->nl = -1;   
  s->wl = 1e-10;

  s->blflag = 0;
  s->blfilter = 1e12;
  s->nfilterflag = 0;
  s->use_rl_filter = 0;

  /*
   * Noise estimation parameters
   */
  s->noise_estimation_flag = 0;
  strcpy( s->qnoisefile, "None" );
  s->gfilter_flag = 0;
  strcpy( s->gfilterfile, "None" );

  /*
   *  Plotting parameters
   */
  s->section = -1;
  s->theta = 0.0;
  s->r = 0.0;
  s->r2 = 0.0;
  
}


/*
 * Write settings to file
 */
void writeConfig(settings * s){

  printf("outpath = %s\n", s->outpath);
  // printf("input correlation function = %s\n", s->correlationfile);
  // printf("input q sampling points = %s\n", s->qfile );
  printf("tag = %s\n", s->tag);
  //  printf("Samples q polar angle (theta_q) = %d\n", s->nthq );
  //  printf("Samples q radial                = %d\n", s->nq );
  printf("Samples r radial                = %d\n", s->nr );
  // printf("qmax                            = %g\n", s->qmax );
  printf("rmax                            = %g\n", s->rmax );
  printf("nl - number of l values         = %d\n", s->nl );
  printf("Section                         = %d\n", s->section );
  printf("Theta                           = %g\n", s->theta );
  printf("r                               = %g\n", s->r );
  printf("r2                              = %g\n", s->r2 );
  //printf("wavelength (metres)             = %g\n", s->wl );
  //printf("noise estimation flag           = %d\n", s->noise_estimation_flag);
  // printf("qnoise file                     = %s\n", s->qnoisefile );
  //printf("blflag                          = %d\n", s->blflag );
  //printf("blfilter                        = %g\n", s->blfilter );
  //printf("nfilterflag                     = %d\n", s->nfilterflag );
  //printf("use_rl_filter                   = %d\n", s->use_rl_filter );
  //printf("input correlation sigma function = %s\n", s->correlation_sigma_file);
  //printf("gfilter_flag                    = %d\n", s->gfilter_flag );
  //printf("gfilter file                     = %s\n", s->gfilterfile );

  FILE * outfile = fopen( s->logname, "a" );

  if (outfile != NULL) {  
    fprintf(outfile,"outpath = %s\n", s->outpath);
    //    fprintf(outfile,"input correlation function = %s\n", s->correlationfile);
    //    fprintf(outfile,"input q sampling points = %s\n", s->qfile );
    fprintf(outfile,"tag = %s\n", s->tag);
    //    fprintf(outfile,"Samples q polar angle (theta_q) = %d\n", s->nthq );
    //  fprintf(outfile,"Samples q radial                = %d\n", s->nq );
    fprintf(outfile,"Samples r radial                = %d\n", s->nr );
    fprintf(outfile,"Section                         = %d\n", s->section );
    fprintf(outfile,"Theta                           = %g\n", s->theta );
    fprintf(outfile,"r                               = %g\n", s->r );
    fprintf(outfile,"r2                              = %g\n", s->r2 );
    //  fprintf(outfile,"qmax                            = %g\n", s->qmax );
    fprintf(outfile,"rmax                            = %g\n", s->rmax );
    fprintf(outfile,"nl - number of l values         = %d\n", s->nl );
    // fprintf(outfile,"wavelength (metres)             = %g\n", s->wl );
    // fprintf(outfile,"noise estimation flag           = %d\n", s->noise_estimation_flag);
    // fprintf(outfile,"qnoise file                     = %s\n", s->qnoisefile );
    // fprintf(outfile,"blflag                          = %d\n", s->blflag );
    // fprintf(outfile,"blfilter                        = %g\n", s->blfilter );
    // fprintf(outfile,"nfilterflag                     = %d\n", s->nfilterflag );
    // fprintf(outfile, "use_rl_filter                   = %d\n", s->use_rl_filter );
    // fprintf(outfile,"input correlation sigma function = %s\n", s->correlation_sigma_file);
    // fprintf(outfile, "gfilter_flag                    = %d\n", s->gfilter_flag );
    // fprintf(outfile, "gfilter file                     = %s\n", s->gfilterfile );
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

  //printf("DEBUG in readConfig outpath = %s\n", s->outpath);

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
  

  if (strcmp(label,"outpath") == 0) {
    strcpy(s->outpath, value);
    if ( (s->outpath[strlen(s->outpath)-1] != '/')
	 || (s->outpath[strlen(s->outpath)-1] != '\\') ){
      strcat( s->outpath, "/");
    }
  }
  else if (strcmp(label,"tag") == 0) {
    strcpy(s->tag, value);
  }
  else if (strcmp(label,"correlationfile") == 0) {
    strcpy(s->correlationfile, value);
    pa = value;
    fp = fopen( pa,   "r");
    if (fp == NULL){
      printf("Correlation file could not be opened - check the directory/name given in the config file\n");
    }
    fclose( fp );
  }
  else if (strcmp(label,"correlation_sigma_file") == 0) {
    strcpy(s->correlation_sigma_file, value);
    pa = value;
    fp = fopen( pa,   "r");
    if (fp == NULL){
      printf("Correlation sigma file could not be opened - check the directory/name given in the config file\n");
    }
    fclose( fp );
  }
  else if (strcmp(label,"qfile") == 0) {
    strcpy(s->qfile, value);
    pa = value;
    fp = fopen( pa,   "r");
    if (fp == NULL){
      printf("qfile could not be opened - check the directory/name given in the config file\n");
    }
    fclose( fp );
  }
  else if (strcmp(label,"qnoisefile") == 0) {
    strcpy(s->qnoisefile, value);
    pa = value;
    fp = fopen( pa,   "r");
    if (fp == NULL){
      printf("qnoisefile could not be opened - check the directory/name given in the config file\n if noise_estimation was not intended, please remove the qnoisefile line from the config file\n ");
    }
    fclose( fp );
  }
  else if (strcmp(label,"gfilterfile") == 0) {
    strcpy(s->gfilterfile, value);
    pa = value;
    fp = fopen( pa,   "r");
    if (fp == NULL){
      printf("gfilterfile could not be opened - check the directory/name given in the config file\n if noise_estimation was not intended, please remove the qnoisefile line from the config file\n ");
    }
    fclose( fp );
  }
  else if (strcmp(label,"nthq") == 0) {
    s->nthq = atoi(value); 
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
  else if (strcmp(label,"wavelength") == 0) {
    s->wl = atof(value); 
  }
  else if (strcmp(label,"rmax") == 0) {
    s->rmax = atof(value); 
  }
  else if (strcmp(label,"nl") == 0) {
    s->nl = atoi(value);  
  }
  else if (strcmp(label,"section") == 0) {
    s->section = atoi(value);  
  }
  else if (strcmp(label,"theta") == 0) {
    s->theta = atof(value); 
  }
  else if (strcmp(label,"r") == 0) {
    s->r = atof(value); 
  }
  else if (strcmp(label,"r2") == 0) {
    s->r2 = atof(value); 
  }
  else if (strcmp(label,"blflag") == 0) {
    s->blflag = atoi(value); 
  }
  else if (strcmp(label,"blfilter") == 0) {
    s->blfilter = atof(value); 
  }
  else if (strcmp(label,"gfilter_flag") == 0) {
    s->gfilter_flag = atoi(value); 
  }
  else if (strcmp(label,"nfilterflag") == 0) {
    s->nfilterflag = atoi(value); 
  }
  else if (strcmp(label,"use_rl_filter") == 0) {
    s->use_rl_filter = atoi(value); 
  }
  else if (strcmp(label,"noise_estimation_flag") == 0) {
    s->noise_estimation_flag = atoi(value); 
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


