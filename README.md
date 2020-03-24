# py2PADF

This Project contains code to the Pair-Angle Distribution Function (PADF) from a set of diffraction patterns. As an intermediate step it calculates the angular cross-correlation function.

WARNING: This is development code. NOT intended for general release. When used correctly it works, but errors in input parameters do not produce meaningful error messages.

The code is made available as a record of the software used to analyse data in the following publications:

Martin et al. "..." Small 2020 (accepted)

For details about the data analysis in this publication, see DATA_ANALYSIS_STEPS file.


**************************************
padf_0.1 - compile & install
**************************************

padf_0.1 calculates the PADF and angular intensity cross-correlation functions from a set of diffraction patterns. It is written with a Python frontend and C backend. The Python code processes input parameters and performs basic processing on the diffraction patterns and correlation functions. The C programs (padf, padfcorr, padfXcorr, padfplot) are separate programs that are run by the Python code (padf_0.1). 

The C programs depend on the Gnu Scientific Library (GSL) and FFTW libraries. No version control issues have been encountered so far with these libraries.

Makefiles for C programs are written with the gcc compiler. 

Compile the C-code from a terminal (linux or OSX).  Change directory to one containing the C-code (padf or padfplot or padfcorr padfXcorr) and type

make clean

make

Repeat for all the C programs (The need to be compiled independently). 


**************************************
padf_0.1  Executing the code
**************************************

From the directory contains the padf_0.1 file:

python padf_0.1 -c <configuration file name>

The configuration file is a text file that contains parameter values.

