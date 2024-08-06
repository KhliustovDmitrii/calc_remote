/****************************\
 * Some useful mathematical *
 * functions                *
 * **************************/
 
void * joint_fdfun_wrapper(void *);

double bessj0( double x )
/*------------------------------------------------------------*/
/* PURPOSE: Evaluate Bessel function of first kind and order  */
/*          0 at input x                                      */
/*------------------------------------------------------------*/
{
   double ax,z;
   double xx,y,ans,ans1,ans2;

   if ((ax=fabs(x)) < 8.0) {
      y=x*x;
      ans1=57568490574.0+y*(-13362590354.0+y*(651619640.7
         +y*(-11214424.18+y*(77392.33017+y*(-184.9052456)))));
      ans2=57568490411.0+y*(1029532985.0+y*(9494680.718
         +y*(59272.64853+y*(267.8532712+y*1.0))));
      ans=ans1/ans2;
   } else {
      z=8.0/ax;
      y=z*z;
      xx=ax-0.785398164;
      ans1=1.0+y*(-0.1098628627e-2+y*(0.2734510407e-4
         +y*(-0.2073370639e-5+y*0.2093887211e-6)));
      ans2 = -0.1562499995e-1+y*(0.1430488765e-3
         +y*(-0.6911147651e-5+y*(0.7621095161e-6
         -y*0.934935152e-7)));
      ans=sqrt(0.636619772/ax)*(cos(xx)*ans1-z*sin(xx)*ans2);
   }
   return ans;
}

double log_lin(double x)
{
    if(x > 1)
       return 1 + log(x);
    else 
    {
       if(x > -1)
          return x;
       else
          return -1 - log(-x);
    }
}

void cube_spline(spline *s,double x0, double x1,double y0,double y1,double dy0,double dy1)
{
     s->a[3] = ( 2*(x0*dy0-y0-x1*dy1+y1)+(dy0-dy1)*(x0+x1) )/((x0-x1)*(7*x0*x0+10*x0*x1+7*x1*x1) );
     s->a[2] = ( dy0-dy1-3*(x0-x1)*(x0+x1)*s->a[3] )/( 2.*(x0-x1) );
     s->a[1] = ( (y0-y1) - (x0*x0*x0-x1*x1*x1)*s->a[3] - (x0-x1)*(x0+x1)*s->a[2] )/( x0-x1 );
     s->a[0] = y0 - x0*x0*x0*s->a[3] - x0*x0*s->a[2] - x0*s->a[1];
}

void cube_spline0(spline *s,double x, double y0,double y1,double dy0,double dy1)
{
     s->a[0] = y0;
     s->a[1] = dy0;
     s->a[2] = ((3*(y1-y0)/x)  - (2*dy0+dy1))/x;
     s->a[3] = ((dy0+dy1) - (2*(y1-y0)/x))/x/x;
}

void compute_gradient(double *x, double *y_ini, double **Jacob, int jdim_1, int jdim_2, int *logpar, 
                      geometry geo, int nlay, int freq_num_td, int freq_num_fd, double *freqs_td, double *freqs_fd, int* tchns,
                      int base_chn, int num_channels, int spec_len, int num_pol_lay, int *ind_pol_lay, double complex *impulse_spec,
                      int *gradpars)
{
   int i, j, k, rc, xsize;
   
   double **x_delta, **y_delta, *div_arr;
   pthread_t threads[NTHREADS];
   thread_arg_struct *thread_args;
   
   x_delta = (double **)malloc(jdim_2*sizeof(double *));
   div_arr = (double *)malloc(jdim_2*sizeof(double));
   thread_args = (thread_arg_struct *)malloc(jdim_2*sizeof(thread_arg_struct));
   for(i = 0; i < jdim_2; i++)
   {
      x_delta[i] = (double *)malloc(jdim_2*sizeof(double));
      for(j = 0; j < jdim_2; j++)
      {
         x_delta[i][j] = x[j];
         if(i==j)
         {
            x_delta[i][j] = (1. - DELTA)*x[j];
            div_arr[i] = log(x_delta[i][j]/x[j]);
            if(logpar[i] == 0)
            {
               x_delta[i][j] = x[j] - 0.1;
               div_arr[i] = 0.1;
            }
         }
      }
   }
   
   y_delta = (double **)malloc(jdim_2*sizeof(double *));
   for(i = 0; i < jdim_2; i++)
   {
      y_delta[i] = (double *)malloc((jdim_1+10)*sizeof(double));
      for(j = 0; j < jdim_1+10; j++) y_delta[i][j] = 0;
   }
   
   for(i = 0; i < jdim_2; i++)
   {
      if(gradpars[i]==1)
      {
         thread_args[i].x = x_delta[i];
         thread_args[i].y = y_delta[i];
         thread_args[i].geo = &geo;
         thread_args[i].nlay = nlay;
         thread_args[i].freq_num_td = freq_num_td;
         thread_args[i].freq_num_fd = freq_num_fd;
         thread_args[i].base_chn = base_chn;
         thread_args[i].num_channels = num_channels;
         thread_args[i].spec_len = spec_len;
         thread_args[i].num_pol_lay = num_pol_lay;
         thread_args[i].freqs_td = freqs_td;
         thread_args[i].freqs_fd = freqs_fd;
         thread_args[i].tchns = tchns;
         thread_args[i].cole_pars = &(x_delta[i][2*nlay]);
         thread_args[i].ind_pol_lay = ind_pol_lay;
         thread_args[i].impulse_spec = impulse_spec;
      }
   }
   if(NTHREADS>1 && jdim_2>1)
   {
      
      j = 0;
      while(j < jdim_2)
      {
             i = 0;
        while(i<NTHREADS)
        { 
           if(gradpars[j]==1)
           {
                
            rc = pthread_create(&threads[i], 0, joint_fdfun_wrapper, (void *) &thread_args[j]);
            i++;
            
           }
           j++;
           if(j == jdim_2) break;
        }

        for (k=0; k<i; k++)
         rc = pthread_join(threads[k], NULL);
      }
   } 
   else
   {
      for(j = 0; j < jdim_2; j++)
      {
           if(gradpars[j]==1)
             joint_fdfun_wrapper((void *) &thread_args[j]); 
      }  
        
   }
   
   for(i = 0; i < jdim_1; i++)
   {
      for(j = 0; j < jdim_2; j++)
      {
         if(gradpars[j]==1)
         {
            switch(DATA_TRANSFORM)
            {
               case 0:
                  Jacob[i][j] = (y_delta[j][i] - y_ini[i])/div_arr[j];
                  break;
               case 1:
                  Jacob[i][j] = log(fabs(y_delta[j][i]/y_ini[i]))/div_arr[j];
                  break;
               case 2:
                  Jacob[i][j] = (log_lin(y_delta[j][i]) - log_lin(y_ini[i]))/div_arr[j];
                  break;
            }
         }
         else
         Jacob[i][j] = 0;
         
         if(isnan(Jacob[i][j])) Jacob[i][j] = 0;
      }
   }
   
   for(i = 0; i < jdim_2; i++) free(x_delta[i]);
   free(x_delta);
   for(i = 0; i < jdim_2; i++) free(y_delta[i]);
   free(y_delta);
   free(div_arr);
   free(thread_args);
}

double primField(double hd,
                 double vd)
{

    double E[10], R[4], RR[10], Hp[4], k, MR, Ampl, M;
    R[1]=hd; R[2]=0; R[3]=vd; //towed cable = 39.11
    //R[1]=33; R[2]=0; R[3]=21; //towed cable = 39.11
    //M[1]=0; M[2]=0; M[3]=1.e9;
    M = 1;
    E[1]=1; E[2]=0; E[3]=0;
    E[4]=0; E[5]=1; E[6]=0;
    E[7]=0; E[8]=0; E[9]=1;
    RR[1]=R[1]*R[1]; RR[2]=R[1]*R[2]; RR[3]=R[1]*R[3];
    RR[4]=R[2]*R[1]; RR[5]=R[2]*R[2]; RR[6]=R[2]*R[3];
    RR[7]=R[3]*R[1]; RR[8]=R[3]*R[2]; RR[9]=R[3]*R[3];
    MR=R[1]*R[1]+R[2]*R[2]+R[3]*R[3];
    k = M / MR / (sqrt(MR)) / 4 / M_PI;
    Hp[1] = (3*RR[3] / MR-E[3]) * k;
    Hp[2] = (3*RR[6] / MR-E[6]) * k;
    Hp[3] = (3*RR[9] / MR-E[9]) * k;

    Ampl = sqrt(Hp[1] * Hp[1] + Hp[2] * Hp[2] + Hp[3] * Hp[3]);
    return Ampl;
}
