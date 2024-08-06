#ifndef _KALMAN_H_

#define _KALMAN_H_

#include <complex.h>

#include "globals.h"

void proc(int, double *, double *, double, double *, double);
int flinversion_fixed(geometry, int, double *, double *, int, int *, double *, double *, double *, int, double complex *,
                int, int *, int, int, int, double *, double *, int *, double *, double *, double *, double *, double *); 
int flinversion_free(geometry, int, double *, double *, int, int *, double *, double *, double *, int, double complex *,
                int, int *, int, int, int, double *, double *, int *, double *, double *, double *, double *, double *);
#endif
