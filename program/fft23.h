#ifndef _FFT23_H_

#define _FFT23_H_

#include <complex.h>

typedef struct
{
    int n;
    int order2;
    int order3;
    int *inv_index;
    double complex *wkn;
    double complex *xn;
    double complex *fn;
} FFT;

FFT *FFT_new(int);
void FFT_free(FFT *);
void FFT_free_full(FFT *);
void fft_pro(FFT *,int);


#endif
