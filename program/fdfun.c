/*******************************\
 * Functions for solving       *
 * forward problem             *
\*******************************/

#include "globals.h"

double bessj0(double);
void cube_spline(spline *,double, double, double, double, double, double);
void cube_spline0(spline *, double, double, double, double, double);

double complex PartSum(double n0,double hh,double complex n1,double r,double complex Imp)
{
    double complex s;
    
    if(fabs(r)>.001)
        s = bessj0(n0*r);
    else 
        s  = 1;

    double complex A = exp(-n0*hh)*(n1-n0*Imp)*n0*n0*.25/(n1+n0*Imp)/M_PI;
    s = A*s;
    return s;
}

double complex Impedance(int n, double n0, double om, double complex *rho, double *dep) 
{
    int i, m;
    double complex ni,nim1;
    double complex Imp;
    double dpth;

    Imp = 1;
    nim1 = csqrt(n0*n0 - I*om*mu0/rho[n-1]);
    dpth = 0;
    m = n - 1;

    for(i=m;i>0;i--)
    {   
       dpth+=dep[i-1];
       ni = nim1;
       nim1 = csqrt(n0*n0 - I*om*mu0/rho[i-1]);
       Imp = ctanh(nim1*dpth+catanh(nim1/ni*Imp));
       dpth = 0;
    }
    return Imp;
}

double complex integral(int n,double hh,double r,double complex *rho,double *dep,double f)
{
    double complex PS;     
    double complex intl = 0;
    double dn0;
    double n0=0;
    double complex n1,c;
    double sigma = 1./rho[0];
    double complex Imp;
    double om = f*2*M_PI;
    c = I*(-om*sigma*mu0);   

    double VAL = .001;
    for(n0=VAL,dn0=VAL;n0<1;n0+=dn0) 
    {
        n1 = csqrt(n0*n0+c);
        Imp = Impedance(n,n0,om,rho,dep);
        PS  = PartSum(n0,hh,n1,r,Imp);
        if(isnan(creal(PS)))
        {
           if(isnan(cimag(PS)))
             PS = 0 + 0*I;
           else
             PS = 0 + cimag(PS)*I;
        } 
        else 
        {
            if(isnan(cimag(PS)))
              PS = creal(PS) + 0*I;
 
        }
        intl += dn0*PS;
    }
    return intl;
}


double complex ImHz(int n, double r,double z,double f,double complex *rho, double *dep) 
{
    return integral(n,z,r,rho,dep,f);
}

//#################################################################

void fd_fdfun_no_diff(geometry geo, int nlay, int bfr,double *x, double *cole_pars, int num_pol_lay, int *ind_pol_lay, double *y, double *freqs) 
{
    int i, j, p;
    int freq_num = bfr;
    double complex refl, alp;
    double complex rho[nlay];
    double da = 0;
    double rho_inf, m, c, tau;
    
    if(bfr==1) freq_num = 2;
    if(nlay>1) da = x[nlay];
    
    for(i=0;i<freq_num;i++) 
    {
        if(num_pol_lay==0)
          for(j = 0; j < nlay; j++) rho[j] = x[j];
        else
        {
            p = 0;
            for(j = 0; j < nlay; j++)
            {
               if(p < num_pol_lay)
               {
                 if(ind_pol_lay[p]==j)
                 {
                     switch(POL_MODEL)
                     {
                     case 0: //Cole-Cole
                        rho_inf = cole_pars[POL_PAR_NUM*p];
                        m = 1 - rho_inf/x[j];
                        tau = cole_pars[POL_PAR_NUM*p + 1];
                        c = cole_pars[POL_PAR_NUM*p + 2];
                        alp = pow(freqs[i]*tau, c);
                        alp = alp * cexp(c*I*5*M_PI/2); //Modified Cole-Cole model
                        //alp = alp * cexp(c*I*M_PI/2); //Vanilla Cole-Cole model
                        rho[j] = x[j]*(1 - m*(1 - 1/alp));
                        break;
                     case 1: //Hyperbolic
                        rho[j] = 1/(cole_pars[POL_PAR_NUM*p]*freqs[i] + cole_pars[POL_PAR_NUM*p + 1]) + x[j];
                        break;
                     }
                     p++;
                 } 
                 else
                   rho[j] = x[j];  
                     
               } 
               else
                 rho[j] = x[j];  
            }   
        }
        refl = ImHz(nlay,geo.hor_dist,2*(geo.alt + da)+geo.ver_dist,freqs[i],rho,&(x[nlay + 1]))*I/geo.prim ;
        y[2*i] = -creal(refl);
        y[2*i+1] = cimag(refl);
        if(i>bfr) break;
    }
}

void fd_fdfun(geometry geo, int nlay, int bfr,double *x, double * cole_pars, int num_pol_lay, int * ind_pol_lay, double *y, double *freqs) 
{
    int i, j, p;
    int freq_num = bfr;
    double complex refl, alp;
    double complex rho[nlay];
    double da = 0;
    double rho_inf, m, c, tau;
    
    if(bfr==1) freq_num = 2;
    if(nlay>1) da = x[nlay];
    
    for(i=0;i<freq_num;i++) 
    {
        if(num_pol_lay==0)
            for(j = 0; j < nlay; j++) rho[j] = x[j];
        else
        {
            p = 0;
            for(j = 0; j < nlay; j++)
            {
               if(p < num_pol_lay)
               {  
                  if(ind_pol_lay[p]==j)
                  {
                     switch(POL_MODEL)
                     {
                     case 0: //Cole-Cole
                        rho_inf = cole_pars[POL_PAR_NUM*p];
                        m = 1 - rho_inf/x[j];
                        tau = cole_pars[POL_PAR_NUM*p + 1];
                        c = cole_pars[POL_PAR_NUM*p + 2];
                        alp = pow(freqs[i]*tau, c);
                        alp = alp * cexp(c*I*5*M_PI/2); //Modified Cole-Cole model
                        //alp = alp * cexp(c*I*M_PI/2); //Vanilla Cole-Cole model 
                        rho[j] = x[j]*(1 - m*(1 - 1/alp));
                        break;
                     case 1: //Hyperbolic
                        rho[j] = 1/(cole_pars[POL_PAR_NUM*p]*freqs[i] + cole_pars[POL_PAR_NUM*p + 1]) + x[j];
                        break;
                     }
                     p++;
                  } 
                  else
                    rho[j] = x[j];  
               }
               else
                 rho[j] = x[j];  
            }   
        }
        refl = ImHz(nlay,geo.hor_dist,2*(geo.alt + da)+geo.ver_dist,freqs[i],rho,&(x[nlay + 1]))*I/geo.prim ;
        y[2*i] = creal(refl);
        y[2*i+1] = cimag(refl);
        if(i>bfr) break;
    }
    
    for(i=0;i<bfr;i++)
    {
      if(i<freq_num-1)
      {
         if(TAKE_ABS==0)
            y[2*i+1] = y[2*(i+1)+1] - y[2*i+1];
         else
            y[2*i+1] = fabs(y[2*(i+1)+1] - y[2*i+1]);
      }
      else
        if(bfr > 1)
           y[2*i+1] = y[2*(i-1)+1];       
    }
     
}
//#################################################################

void td_fdfun(geometry geo, int nlay, double *x, double * cole_pars, int num_pol_lay, 
              int * ind_pol_lay, double *y, double *freqs, int freq_num, double complex *impulse_spec, int spec_len) 
{
    int sfr, j;
    spline sr, si;
    
    double *fdfun_spec;
    
    fdfun_spec = (double *) malloc(2*freq_num*sizeof(double));
    
    fd_fdfun_no_diff(geo, nlay, freq_num, x, cole_pars, num_pol_lay, ind_pol_lay, fdfun_spec, freqs);
    
    // clearing high frequencies
    double complex y0 =  fdfun_spec[2*freq_num-1]*I + fdfun_spec[2*freq_num-2];
    double dy0 =  0;
  
    double k = spec_len-1;
    double y1 = 0;
    double dy1 = 0;
    double lj;

    double xx = log(k)-log(freq_num-1);
    cube_spline0(&sr   ,xx ,creal(y0)   ,creal(y1)   ,creal(dy0)   ,creal(dy1)  );
    cube_spline0(&si   ,xx ,cimag(y0)   ,cimag(y1)   ,cimag(dy0)   ,cimag(dy1)  );
    
    
    pthread_mutex_lock(&fft_mutex);
    memset(&ifft,0,sizeof(ifft));
    init_fft(&ifft,spec_len);
    
    for(sfr = 1; sfr < 2*freq_num; sfr+=2)
       ifft.xn[sfr] = impulse_spec[sfr]*(freqs[(sfr-1)/2]*(-fdfun_spec[sfr-1] + I*fdfun_spec[sfr]))*2*M_PI*100;
    
    for(j=2*freq_num;j<=k;j+=2) 
    {
       lj = log(j)-log(2*freq_num);
       ifft.xn[j] = sr.a[3]*lj*lj*lj+sr.a[2]*lj*lj+sr.a[1]*lj+sr.a[0] +
                            I * (si.a[3]*lj*lj*lj+si.a[2]*lj*lj+si.a[1]*lj+si.a[0]);                     
    }
       
    fft_pro(&ifft,1);
    for(sfr = 0; sfr < spec_len; sfr++)
       y[sfr] = creal(ifft.fn[sfr]);

    FFT_free_full(&ifft);
    pthread_mutex_unlock(&fft_mutex);
    free(fdfun_spec);
}

void cut_into_channels(int *tchns, int base_chn, int num_channels, double *y, double *output)
{
   
   int i, j;
   
   for(i = 0; i < num_channels; i++)
   {
      output[i] = 0;
      for(j = tchns[i]; j < tchns[i+1]; j++)
         output[i] += y[base_chn + j];
            
      output[i] = output[i]/(tchns[i+1] - tchns[i]);
   }
}

void chan_td_fdfun(geometry geo, int nlay, double *x, double * cole_pars, int num_pol_lay, int * ind_pol_lay, double *y_chan, double *y_full, double *freqs,
                   int freq_num, double complex *impulse_spec, int spec_len,
                   int* tchns, int base_chn, int num_channels)
{
   td_fdfun(geo, nlay, x, cole_pars, num_pol_lay, ind_pol_lay,y_full, freqs, freq_num, impulse_spec, spec_len);
   cut_into_channels(tchns, base_chn, num_channels, y_full, y_chan);
}

void joint_fdfun(geometry geo, int nlay, double *x, double * cole_pars, int num_pol_lay, int * ind_pol_lay, double *y_full, double *freqs_td,
                   int freq_num_td, double complex *impulse_spec, int spec_len,
                   int* tchns, int base_chn, int num_channels,
                   int freq_num_fd, double *freqs_fd,
                   double *y_joint)
{
   if(num_channels>0)
   {
      td_fdfun(geo, nlay, x, cole_pars, num_pol_lay, ind_pol_lay,y_full, freqs_td, freq_num_td, impulse_spec, spec_len);
      cut_into_channels(tchns, base_chn, num_channels, y_full, &(y_joint[2*freq_num_fd]));
   }
   fd_fdfun(geo, nlay, freq_num_fd, x, cole_pars, num_pol_lay, ind_pol_lay,y_joint, freqs_fd);
}

void * joint_fdfun_wrapper(void * arg)
{
   thread_arg_struct * str = (thread_arg_struct *) arg;
   double *y_full;
   y_full = (double *)malloc(str->spec_len * sizeof(double));
   joint_fdfun(*(str->geo), str->nlay, str->x, str->cole_pars, str->num_pol_lay, str->ind_pol_lay, y_full, str->freqs_td,
                   str->freq_num_td, str->impulse_spec, str->spec_len,
                   str->tchns, str->base_chn, str->num_channels,
                   str->freq_num_fd, str->freqs_fd,
                   str->y);
    free(y_full);
}
