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
  char input_dp_name[1024];
  char input_dp_name2[1024];
  char logname[1024];
  char outpath[1024];
  int output_xsections;
  int Xflag;

  /*
   *  Diffraction parameters
   */
  int cx;  // x position of centre of diffraction pattern
  int cy;  // y position of centre of diffraction pattern

  /*
   *  Polar coordinate parameters
   */
  int nth;  // number of polar angle samples
  int nq;   // number of radial q samples    //NOT USED
  int nr;   // number of real space samples
  double qmax;    // maximum q in interpolation
  double rmax;    // maximum r to include in the calculation
  
  /*
   *  Scale parameters
   */
  double  pixel_width; // width of a pixel in metres (used to define output array size)
  double  detector_z; //  detector distance in metres
  double  wavelength; // in metres

  /*
   *  array dimensions (if not set by config, then set in code)
   */
  int nx;
  int ny;
  double pwx;
  double pwy;

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

void default_settings( settings * s, int nx );
#endif
