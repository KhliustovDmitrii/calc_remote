#ifndef _MU_H_

#define _MU_H_

#include <complex.h>

#include "globals.h"

double bessj0(double);
double log_lin(double);
void cube_spline(spline *,double, double, double, double, double, double);
void cube_spline0(spline *, double, double, double, double, double);
void compute_gradient(double *, double *, double **, int, int, int *, geometry, int, int, int, double *, double *, int *, int, 
                      int, int, int, int *, double complex *, int *);
double primField(double, double);

#endif
