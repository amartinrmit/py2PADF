/*
 * settings.h : part of padf
 *
 * A.V. Martin March 2015
 */


#ifndef SETTINGS_H
#define SETTINGS_H


#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <ctype.h>

typedef struct{
    

  /*
   *  name of input file; the output log file
   */
  char tag[1024];
  char correlationfile[1024];
  char correlation_sigma_file[1024];
  char qfile[1024];
  char logname[1024];
  char outpath[1024];
  int output_padf;

  /*
   *  Polar coordinate parameters
   */
  int nthq;  // number of polar angle samples
  int nq;   // number of radial q samples
  int nr;   // number of real space samples
  double qmax;    // maximum q in interpolation
  double rmax;    // maximum r to include in the calculation
  double wl; // wavelength in metres
  
  /*
   *  Spherical harmonic parameters
   */
  int nl;
  int nlmin;

  /*
   *  B_l filter parameters
   */
  int blflag;
  double blfilter;

  // correlation sigma flag
  int nfilterflag;

  // r_l_filter variables
  int use_rl_filter;

  /*
   *  parameters for noise estimation
   */
  int noise_estimation_flag;
  char qnoisefile[1024];
  
  // gaussian width
  int gfilter_flag;
  char gfilterfile[1024];


 } settings;


/*
 *  input file names, experiment parameters, geometery
 *  currently hard coded into settings.c
 */
void initialize_settings( settings * s);

/*
 * Write settings to file
 */
void writeConfig(settings * s);

/*
 *  Read a config file
 */
void readConfig(char * fname, settings * s);
void parseLine(char * line, settings * s);
void breakLine(char * line, char * label, char * value);
void convertToLower(char * str);

/*
 *  read the name of configuration file from command line arguments
 */
void parse_config_name( int argc, char * argv[], char * cname );

char *trimwhitespace(char *str);
#endif
