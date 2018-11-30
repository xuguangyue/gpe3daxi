/**
 *    BEC-GP-OMP codes are developed and (c)opyright-ed by:
 *
 *    Luis E. Young-S., Sadhan K. Adhikari
 *    (UNESP - Sao Paulo State University, Brazil)
 *
 *    Paulsamy Muruganandam
 *    (Bharathidasan University, Tamil Nadu, India)
 *
 *    Dusan Vudragovic, Antun Balaz
 *    (Scientific Computing Laboratory, Institute of Physics Belgrade, Serbia)
 *
 *    Public use and modification of this code are allowed provided that the
 *    following three papers are cited:
 *
 *    [1] L. E. Young-S. et al., Comput. Phys. Commun. 204 (2016) 209.
 *    [2] P. Muruganandam, S. K. Adhikari, Comput. Phys. Commun. 180 (2009) 1888.
 *    [3] D. Vudragovic et al., Comput. Phys. Commun. 183 (2012) 2021.
 *
 *    The authors would be grateful for all information and/or comments
 *    regarding the use of the code.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include <omp.h>
#include <time.h>

#define MAX(a, b) (a > b) ? a : b
#define RMS_ARRAY_SIZE     4
#define MAX_FILENAME_SIZE  256

#define BOHR_RADIUS         5.2917720859e-11

char *output, *initout, *rmsout, *Nstpout, *Npasout, *Nrunout, *dynaout, *tempout;
long outstprho, outstpz, outstpt, outstpwf;

int opt;
long Na;
long Nrho, Nz;
long Nrho2, Nz2;
long Nstp, Npas, Nrun;
double drho, dz;
double drho2, dz2;
double dt;
double G0, Gpar, G;
double aho, as;
double vnu, vlambda;
double vnu2, vlambda2;
double ampr, ampz, zamp, freq;
double par;
double pi;

double *rho, *z;
double *rho2, *z2;
double **pot;

double complex Arho0, Az0, Arho0r, Az0r, Arho, Az, dArho;
double complex *calpharho, *calphaz;
double complex *cgammarho, *cgammaz;

void readpar(void);
void init(double complex **, double **);
double randn(double, double);
void gencoef(void);
void calcnorm(double *, double complex **, double **, double **);
void calcmuen(double *, double *, double complex **, double **, double **, double **, double **, double **, double **, double **, double **, double **, double **);
void calcrms(double *, double complex **, double **, double **, double **, double **);
void calcnu(double complex **);
void calclurho(double complex **, double complex **);
void calcluz(double complex **, double complex **);
void outpsirhoz(double complex **, FILE *);
void outdenrho(double complex **, double *, FILE *);
void outdenz(double complex **, double *, FILE *);

extern double simpint(double, double *, long);
extern void diff(double, double *, double *, long);

extern double *alloc_double_vector(long);
extern double complex *alloc_complex_vector(long);
extern double **alloc_double_matrix(long, long);
extern double complex **alloc_complex_matrix(long, long);
extern void free_double_vector(double *);
extern void free_complex_vector(double complex *);
extern void free_double_matrix(double **);
extern void free_complex_matrix(double complex **);

extern int cfg_init(char *);
extern char *cfg_read(char *);
