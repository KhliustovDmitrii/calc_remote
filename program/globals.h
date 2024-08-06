#ifndef _GLOBAL_H_

#define _GLOBAL_H_

#include <complex.h>

typedef struct {
     double a[4];
} spline;

typedef struct {
    double coord;
    double coord_y;
    double hor_dist;
    double ver_dist;
    double relief;
    double alt;
    double prim;
    double *w;
} geometry;

typedef struct{
   double *x;
   double *y;
   geometry *geo;
   int nlay;
   int freq_num_td;
   int freq_num_fd;
   int base_chn;
   int num_channels;
   int spec_len;
   int num_pol_lay;
   double *freqs_td;
   double *freqs_fd;
   int* tchns;
   double *cole_pars;
   int *ind_pol_lay;
   double complex *impulse_spec;
} thread_arg_struct;

typedef struct{
  int freq_num_fd;
  double freqs_fd[256];
  int lay_num;
  double first_thick;
  double step;
  double ERR_INI;
  double COR_INI;
  int first_mes_pos;
  int last_mes_pos;
  int first_chan_pos;
  int last_chan_pos;
  int hor_dist_pos;
  int ver_dist_pos;
  int alt_pos;
  int time_size;
  int pol_num;
  int pol_inds[64];
  double base_freq;
  int freq_num_td;
  int spec_len;
  int base_chn;
  int num_channels;
  int tchns[64];
  double freqs_td[256];
} config_struct;

#define mu0 (M_PI*.4e-6)
#define DELTA   .1

#define MIN_RES 0.01
#define MAX_RES 20000.
#define RES_INI 1500
#define AVERAGE 4
#define MAX_ITER 3
#define STOP_VAL 1.0

#define TD_ONLY 0

#define INITIAL_METHOD 0   //0 - RES_INI half space, 1 - Half space, 2 - Noise adaptive half space, 3 - KDE

#define POL_MODEL 0 //0 - Cole-Cole, 1 - Hyperbolic 1/(a omega + b) + c
#define POL_PAR_NUM 3 // 3 for Cole-Cole, 2 for Hyperbolic

#define FD_COMP_ORDER 0 //0 - first real, second imag, 1 - vice versa
#define TAKE_ABS 1 //Whether to take fabs of measurements and differences. 0 for no, 1 for yes

#define DATA_TRANSFORM 1 //0 - linear, 1 - logarithmic, 2 - piecewise log-lin
#define PROG_REGIME 0 //0 - local, 1 - server execution

#define LAYERWISE 0 //0 - regular run, 1 - script run for polarization in each layer

#define BUF_SIZE 4096

#define MAX_FREE_LAYERS 6

#define GREEDY 0
#define VERBOSE 1

double *upper;
double *lower;

int NTHREADS;

FFT fft,ifft,fftc,ifftd;

geometry geo, geo1, geo_ini;
config_struct configs;

static pthread_mutex_t fft_mutex = PTHREAD_MUTEX_INITIALIZER;

#endif
