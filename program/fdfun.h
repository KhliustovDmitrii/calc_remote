#ifndef _FDFUN_H_

#define _FDFUN_H_

#include <complex.h>

#include "globals.h"

double complex PartSum(double, double, double complex, double, double complex);
double complex Impedance(int, double, double, double complex, double);
double complex integral(int, double, double, double complex, double, double);
double complex ImHz(int, double, double, double, double complex, double);
void fd_fdfun_no_diff(geometry, int, int, double *, double *, int, int *, double *, double *);
void fd_fdfun(geometry, int, int, double *, double *, int, int *, double *, double *);
void td_fdfun(geometry, int, double *, double *, int, int *, double *, double *, int, double complex *, int);
void cut_into_channels(int *, int, int, double *, double *);
void chan_td_fdfun(geometry, int, double *, double *, int, int *, double *, double *, double *, int, double complex *, int, int *, int, int);
void joint_fdfun(geometry, int, double *, double *, int, int *, double *, double *, int, double complex *, int, int *, int, int, int, double *, double *);
void * joint_fdfun_wrapper(void *);

#endif
