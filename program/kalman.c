/**********************************\
 * Functions for Kalman inversion *
\**********************************/

void joint_fdfun(geometry, int, double *, double *, int, int *, double *, double *, int, double complex *, int, int *, int, int, int, double *, double *);
void compute_gradient(double *, double *, double **, int, int, int *, geometry, int, int, int, double *, double *, int *, int, 
                      int, int, int, int *, double complex *, int *);


void proc(int n, double *x, double *S, double z, double *h, double sg2) 
{
    int i,j,k;
    double d[2],bk,ck,dz;
    
    double *f, *e;
    f = (double *) malloc((n + 10)*sizeof(double));
    e = (double *) malloc((n + 10)*sizeof(double));
    
    for(i = 0; i < n+10; i++)
    {
       f[i] = 0;
       e[i] = 0;
    }
    
    d[0] = sg2;
    for(i=0;i<n;i++) 
    {
        f[i] = 0;
        for(j=0;j<n;j++)
            f[i]+=S[j*n+i]*h[j];
    }
    memset(e,0,sizeof(e));
    for(k=0;k<n+1;k++) 
    {
        d[1] = d[0] + f[k]*f[k];
        bk = sqrt(d[0]/d[1]);
        ck = f[k]/sqrt(d[0]*d[1]);
        for(i=0;i<n;i++) 
        {
            double tmp = S[i*n+k]*f[k];
            S[i*n+k] = bk*S[i*n+k]-ck*e[i];
            e[i] += tmp;
        }
        d[0] = d[1];
    }
    dz = z;
    for(i=0;i<n;i++) dz -= h[i]*x[i];
    dz/=d[0];
    for(i=0;i<n;i++) x[i] += e[i]*dz;
    
    free(f);
    free(e);
}

int flinversion_fixed(geometry geo,
                int nlay,
                double *x_ini,
                double *cole_ini,
                int num_pol_lay,
                int *ind_pol_lay,
                double *dpth,
                double *y_ini,
                double *y_ini_full,
                int freq_num_td, 
                double complex *impulse_spec,
                int spec_len,
                int* tchns,
                int base_chn,
                int num_channels,
                int freq_num_fd,
                double *y_mes,
                double *residual,
                int *up,
                double *S,
                double *freqs_fd,
      double *freqs_td,
      double *upper,
      double *lower) 
{

    int lay_num = nlay;
    double y1[num_channels+2*freq_num_fd];
    double dx[2*lay_num + 1],x0[2*lay_num + 1],x1[2*lay_num + 1],xini[2*lay_num + 1];
    double dcole[POL_PAR_NUM*num_pol_lay],cole0[POL_PAR_NUM*num_pol_lay],cole1[POL_PAR_NUM*num_pol_lay],coleini[POL_PAR_NUM*num_pol_lay];
    double dfull[2*lay_num + 1 + POL_PAR_NUM*num_pol_lay], fullini[2*lay_num + 1 + POL_PAR_NUM*num_pol_lay];
    int logpar[2*lay_num + 1 + POL_PAR_NUM*num_pol_lay], gradpars[2*lay_num + 1 + POL_PAR_NUM*num_pol_lay];
    int charge = 1;
    
    double Jacob_red[num_channels+2*freq_num_fd][lay_num + 1 + POL_PAR_NUM*num_pol_lay];

    double res = 0;
    double k = 10;
    double div;
    int i,j;
    int jdim_1, jdim_2;
    
    
    
    double **Jacob;
    Jacob = (double **) malloc((num_channels+2*freq_num_fd+10)*sizeof(double *));
    for(i = 0; i < num_channels+2*freq_num_fd+10; i++)
    {
       Jacob[i] = (double *) malloc((2*lay_num + 10 + POL_PAR_NUM*num_pol_lay)*sizeof(double));
       for(j = 0; j < 2*lay_num + 10 + POL_PAR_NUM*num_pol_lay; j++) Jacob[i][j] = 0;
    }

    memset(dx,0,sizeof(dx));
    memset(dfull,0,sizeof(dfull));
    for(i=0;i<nlay + 1;i++) 
    {
        xini[i] = x_ini[i];
        fullini[i] = xini[i];
        if(i<nlay-1)
            fullini[i+nlay+1] = xini[i+nlay + 1] = x0[i+nlay + 1] = x1[i+nlay + 1] = dpth[i];
    }
    
    for(i=0;i<POL_PAR_NUM*num_pol_lay;i++)
    {
      coleini[i] = cole_ini[i];
      fullini[i+2*nlay] = coleini[i];
    }
   
   for(i = 0; i < 2*lay_num + 1 + POL_PAR_NUM*num_pol_lay; i++)
   {
       logpar[i] = 1;
       gradpars[i] = 1;
   }
   
   for(i = lay_num+1; i < 2*lay_num; i++)
       gradpars[i] = 0;
       
   logpar[lay_num] = 0;
   if(POL_MODEL==1) 
       for(i = 2*lay_num + 1; i < 2*lay_num + 1 + POL_PAR_NUM*num_pol_lay; i++)
            logpar[i] = 0;

    if(TAKE_ABS==1)
      for(i=0;i<2*freq_num_fd;i++)
         if(y_mes[i]<.001) y_mes[i] = .001;  

    // first forward calculation for the model
    if(*residual<0) 
    {
        joint_fdfun(geo, nlay, xini, coleini, num_pol_lay, ind_pol_lay, y_ini_full, freqs_td,
                   freq_num_td, impulse_spec, spec_len,
                   tchns, base_chn, num_channels,
                   freq_num_fd, freqs_fd, y_ini);
        *residual = 0;
        for(j=0;j<num_channels+2*freq_num_fd;j++) 
        {
            double val = (y_ini[j]-y_mes[j])*(y_ini[j]-y_mes[j])*geo.w[j]*geo.w[j]/(num_channels+2*freq_num_fd);
            *residual += val;
        }
    }

    // a good enough first approach
    if(*residual<0.01)
        return 0;
        
    jdim_1 = num_channels+2*freq_num_fd;
    if(freq_num_fd==1)
       jdim_2 = 1 + POL_PAR_NUM*num_pol_lay;
    else
       jdim_2 = 2*lay_num + POL_PAR_NUM*num_pol_lay;
   
    // Jacobian matrix calculation
   compute_gradient(fullini, y_ini, Jacob, jdim_1, jdim_2, logpar, 
                      geo, nlay, freq_num_td, freq_num_fd, freqs_td, freqs_fd, tchns,
                      base_chn, num_channels, spec_len, num_pol_lay, ind_pol_lay, impulse_spec, gradpars);
                      
   if(VERBOSE>1)
   {
      printf("\n-----------------------------H-------------------------------:\n");
      for(i=0; i<jdim_1; i++)
      {
         for(j = 0; j < jdim_2; j++)
            printf("%f ",Jacob[i][j]);
         printf("\n");
      }
      printf("\n--------------------------------------------------------------\n");            
   }
    
   for(i = 0; i < lay_num + 1; i++)
      for(j = 0; j < num_channels+2*freq_num_fd; j++)
        Jacob_red[j][i] = Jacob[j][i];   
    
   for(i = 0; i < POL_PAR_NUM*num_pol_lay; i++)
      for(j = 0; j < num_channels+2*freq_num_fd; j++) 
        Jacob_red[j][i+lay_num+1] = Jacob[j][i+2*lay_num]; 
    
   switch(DATA_TRANSFORM)
   {
   case 0:
      for(j=0;j<num_channels+2*freq_num_fd;j++)
         proc(nlay+1+POL_PAR_NUM*num_pol_lay,dfull,S, y_ini[j] - y_mes[j],Jacob_red[j],1./(geo.w[j]*geo.w[j]));
      break;
   case 1:
      for(j=0;j<num_channels+2*freq_num_fd;j++)
         proc(nlay+1+POL_PAR_NUM*num_pol_lay,dfull,S, log(fabs(y_ini[j]/y_mes[j])),Jacob_red[j],1./(geo.w[j]*geo.w[j]*y_mes[j]*y_mes[j]));
      break;
   case 2:
      for(j=0;j<num_channels+2*freq_num_fd;j++)
         proc(nlay+1+POL_PAR_NUM*num_pol_lay,dfull,S, log_lin(y_ini[j]) - log_lin(y_mes[j]),Jacob_red[j],1./(geo.w[j]*geo.w[j]));
      break;
    } 

   for(i=0;i<nlay + 1 + POL_PAR_NUM*num_pol_lay;i++) 
   {
        if(i < nlay) x0[i] = xini[i]*exp(-dfull[i]);
        if(i == nlay) x0[i] = xini[i] + dfull[i];
        if(i > nlay) cole0[i - nlay - 1] = coleini[i - nlay - 1]*exp(-dfull[i]);
        if(i <= nlay)
        {
            if(isnan(x0[i])) x0[i] = xini[i];
            if(x0[i]>upper[i]) x0[i] = upper[i];
            if(x0[i]<lower[i]) x0[i] = lower[i];
        } 
        else 
        {
            if(isnan(cole0[i - nlay - 1])) cole0[i - nlay - 1] = coleini[i - nlay - 1];
            if(cole0[i - nlay - 1]>upper[i]) cole0[i - nlay - 1] = upper[i];
            if(cole0[i - nlay - 1]<lower[i]) cole0[i - nlay - 1] = lower[i];
        }
    }

    joint_fdfun(geo, nlay, x0, cole0, num_pol_lay, ind_pol_lay, y_ini_full, freqs_td,
                   freq_num_td, impulse_spec, spec_len,
                   tchns, base_chn, num_channels,
                   freq_num_fd, freqs_fd, y_ini);
                   
    for(j=0;j<num_channels+2*freq_num_fd;j++)
    {
        double val = (y_ini[j]-y_mes[j])*(y_ini[j]-y_mes[j])*geo.w[j]*geo.w[j]/(num_channels+2*freq_num_fd);
        res += val;
    }

    int cntr = 0;
    while(res>*residual*1.01) 
    {
        if(cntr++ > 3) break;
        for(i=0;i<nlay + 1 + POL_PAR_NUM*num_pol_lay;i++) 
        {
          dfull[i]*=.5;
          if(i < nlay) x0[i] = xini[i]*exp(-dfull[i]);
          if(i == nlay) x0[i] = xini[i] + dfull[i];
          if(i > nlay) cole0[i - nlay - 1] = coleini[i - nlay - 1]*exp(-dfull[i]);
          if(i <= nlay)
          {
             if(isnan(x0[i])) x0[i] = xini[i];
             if(x0[i]>upper[i]) x0[i] = upper[i];
             if(x0[i]<lower[i]) x0[i] = lower[i];
          } 
          else 
          {
             if(isnan(cole0[i - nlay - 1])) cole0[i - nlay - 1] = coleini[i - nlay - 1];
             if(cole0[i - nlay - 1]>upper[i]) cole0[i - nlay - 1] = upper[i];
             if(cole0[i - nlay - 1]<lower[i]) cole0[i - nlay - 1] = lower[i];
          }
       }
       res = 0;
       joint_fdfun(geo, nlay, x0, cole0, num_pol_lay, ind_pol_lay, y_ini_full, freqs_td,
                   freq_num_td, impulse_spec, spec_len,
                   tchns, base_chn, num_channels,
                   freq_num_fd, freqs_fd, y_ini);
       for(j=0;j<num_channels+2*freq_num_fd;j++)
       {
            double val = (y_ini[j]-y_mes[j])*(y_ini[j]-y_mes[j])*geo.w[j]*geo.w[j]/(num_channels+2*freq_num_fd);
            res += val;
       }
    }

    if(res>*residual) up[0]++;
    else up[0] = 0;
    *residual = res;
    for(i=0;i<nlay + 1;i++)
        x_ini[i] = x0[i];
    for(i=0;i<POL_PAR_NUM*num_pol_lay;i++)
        cole_ini[i] = cole0[i];
        
    for(i = 0; i < num_channels+2*freq_num_fd+10; i++)
       free(Jacob[i]);

    free(Jacob);
    
    return 1;
}

int flinversion_free(geometry geo,
                int nlay,
                double *x_ini,
                double *cole_ini,
                int num_pol_lay,
                int *ind_pol_lay,
                double *dpth,
                double *y_ini,
                double *y_ini_full,
                int freq_num_td, 
                double complex *impulse_spec,
                int spec_len,
                int* tchns,
                int base_chn,
                int num_channels,
                int freq_num_fd,
                double *y_mes,
                double *residual,
                int *up,
                double *S,
                double *freqs_fd,
      double *freqs_td,
      double *upper,
      double *lower) 
{

    int lay_num = nlay;
    double y1[num_channels+2*freq_num_fd];
    double dx[2*lay_num + 1],x0[2*lay_num + 1],x1[2*lay_num + 1],xini[2*lay_num + 1];
    double dcole[POL_PAR_NUM*num_pol_lay],cole0[POL_PAR_NUM*num_pol_lay],cole1[POL_PAR_NUM*num_pol_lay],coleini[POL_PAR_NUM*num_pol_lay];
    double dfull[2*lay_num + 1 + POL_PAR_NUM*num_pol_lay], fullini[2*lay_num + 1 + POL_PAR_NUM*num_pol_lay];
    int logpar[2*lay_num + 1 + POL_PAR_NUM*num_pol_lay], gradpars[2*lay_num + 1 + POL_PAR_NUM*num_pol_lay];
    int charge = 1;
    
    double Jacob_red[num_channels+2*freq_num_fd][2*lay_num + POL_PAR_NUM*num_pol_lay];

    double res = 0;
    double k = 10;
    double div;
    int i,j;
    int jdim_1, jdim_2;
    
    
    
    double **Jacob;
    Jacob = (double **) malloc((num_channels+2*freq_num_fd+10)*sizeof(double *));
    for(i = 0; i < num_channels+2*freq_num_fd+10; i++)
    {
       Jacob[i] = (double *) malloc((2*lay_num + 10 + POL_PAR_NUM*num_pol_lay)*sizeof(double));
       for(j = 0; j < 2*lay_num + 10 + POL_PAR_NUM*num_pol_lay; j++) Jacob[i][j] = 0;
    }

    memset(dx,0,sizeof(dx));
    memset(dfull,0,sizeof(dfull));
    for(i=0;i<nlay + 1;i++) 
    {
        xini[i] = x_ini[i];
        fullini[i] = xini[i];
        if(i<nlay-1)
            fullini[i+nlay+1] = xini[i+nlay + 1] = x0[i+nlay + 1] = x1[i+nlay + 1] = dpth[i];
    }
    
    for(i=0;i<POL_PAR_NUM*num_pol_lay;i++)
    {
      coleini[i] = cole_ini[i];
      fullini[i+2*nlay] = coleini[i];
    }
   
   for(i = 0; i < 2*lay_num + 1 + POL_PAR_NUM*num_pol_lay; i++)
   {
       logpar[i] = 1;
       gradpars[i] = 1;
   }
       
   logpar[lay_num] = 0;
   if(POL_MODEL==1) 
       for(i = 2*lay_num + 1; i < 2*lay_num + 1 + POL_PAR_NUM*num_pol_lay; i++)
            logpar[i] = 0;

    if(TAKE_ABS==1)
      for(i=0;i<2*freq_num_fd;i++)
         if(y_mes[i]<.001) y_mes[i] = .001;  

    // first forward calculation for the model
    if(*residual<0) 
    {
        joint_fdfun(geo, nlay, xini, coleini, num_pol_lay, ind_pol_lay, y_ini_full, freqs_td,
                   freq_num_td, impulse_spec, spec_len,
                   tchns, base_chn, num_channels,
                   freq_num_fd, freqs_fd, y_ini);
        *residual = 0;
        for(j=0;j<num_channels+2*freq_num_fd;j++) 
        {
            double val = (y_ini[j]-y_mes[j])*(y_ini[j]-y_mes[j])*geo.w[j]*geo.w[j]/(num_channels+2*freq_num_fd);
            *residual += val;
        }
    }

    // a good enough first approach
    if(*residual<0.01)
        return 0;
        
    jdim_1 = num_channels+2*freq_num_fd;
    if(freq_num_fd==1)
       jdim_2 = 1 + POL_PAR_NUM*num_pol_lay;
    else
       jdim_2 = 2*lay_num + POL_PAR_NUM*num_pol_lay;
   
    // Jacobian matrix calculation
   compute_gradient(fullini, y_ini, Jacob, jdim_1, jdim_2, logpar, 
                      geo, nlay, freq_num_td, freq_num_fd, freqs_td, freqs_fd, tchns,
                      base_chn, num_channels, spec_len, num_pol_lay, ind_pol_lay, impulse_spec, gradpars);
    
   if(VERBOSE>11)
   {
      printf("\n-----------------------------H-------------------------------:\n");
      for(i=0; i<jdim_1; i++)
      {
         for(j = 0; j < jdim_2; j++)
            printf("%f ",Jacob[i][j]);
         printf("\n");
      }
      printf("\n--------------------------------------------------------------\n");            
   }
   
    
   switch(DATA_TRANSFORM)
   {
   case 0:
      for(j=0;j<num_channels+2*freq_num_fd;j++)
         proc(2*nlay+POL_PAR_NUM*num_pol_lay,dfull,S, y_ini[j] - y_mes[j],Jacob[j],1./(geo.w[j]*geo.w[j]));
      break;
   case 1:
      for(j=0;j<num_channels+2*freq_num_fd;j++)
         proc(2*nlay+POL_PAR_NUM*num_pol_lay,dfull,S, log(fabs(y_ini[j]/y_mes[j])),Jacob[j],1./(geo.w[j]*geo.w[j]*y_mes[j]*y_mes[j]));
      break;
   case 2:
      for(j=0;j<num_channels+2*freq_num_fd;j++)
         proc(2*nlay+POL_PAR_NUM*num_pol_lay,dfull,S, log_lin(y_ini[j]) - log_lin(y_mes[j]),Jacob[j],1./(geo.w[j]*geo.w[j]));
      break;
    } 

   for(i=0;i<2*nlay + POL_PAR_NUM*num_pol_lay;i++) 
   {
        if(i < 2*nlay) x0[i] = xini[i]*exp(-dfull[i]);
        if(i == nlay) x0[i] = xini[i] + dfull[i];
        if(i >= 2*nlay) cole0[i - 2*nlay] = coleini[i - 2*nlay]*exp(-dfull[i]);
        if(i < 2*nlay)
        {
            if(isnan(x0[i])) x0[i] = xini[i];
            if(x0[i]>upper[i]) x0[i] = upper[i];
            if(x0[i]<lower[i]) x0[i] = lower[i];
        } 
        else 
        {
            if(isnan(cole0[i - 2*nlay])) cole0[i - 2*nlay] = coleini[i - 2*nlay];
            if(cole0[i - 2*nlay]>upper[i]) cole0[i - 2*nlay] = upper[i];
            if(cole0[i - 2*nlay]<lower[i]) cole0[i - 2*nlay] = lower[i];
        }
    }

    joint_fdfun(geo, nlay, x0, cole0, num_pol_lay, ind_pol_lay, y_ini_full, freqs_td,
                   freq_num_td, impulse_spec, spec_len,
                   tchns, base_chn, num_channels,
                   freq_num_fd, freqs_fd, y_ini);
                   
    for(j=0;j<num_channels+2*freq_num_fd;j++)
    {
        double val = (y_ini[j]-y_mes[j])*(y_ini[j]-y_mes[j])*geo.w[j]*geo.w[j]/(num_channels+2*freq_num_fd);
        res += val;
    }

    int cntr = 0;
    while(res>*residual*1.01) 
    {
        if(cntr++ > 3) break;
        for(i=0;i<2*nlay + POL_PAR_NUM*num_pol_lay;i++) 
        {
          dfull[i]*=.5;
          if(i < 2*nlay) x0[i] = xini[i]*exp(-dfull[i]);
          if(i == nlay) x0[i] = xini[i] + dfull[i];
          if(i >= 2*nlay) cole0[i - 2*nlay] = coleini[i - 2*nlay]*exp(-dfull[i]);
          if(i < 2*nlay)
          {
             if(isnan(x0[i])) x0[i] = xini[i];
             if(x0[i]>upper[i]) x0[i] = upper[i];
             if(x0[i]<lower[i]) x0[i] = lower[i];
          } 
          else 
          {
             if(isnan(cole0[i - 2*nlay])) cole0[i - 2*nlay] = coleini[i - 2*nlay];
             if(cole0[i - 2*nlay]>upper[i]) cole0[i - 2*nlay] = upper[i];
             if(cole0[i - 2*nlay]<lower[i]) cole0[i - 2*nlay] = lower[i];
          }
       }
       res = 0;
       joint_fdfun(geo, nlay, x0, cole0, num_pol_lay, ind_pol_lay, y_ini_full, freqs_td,
                   freq_num_td, impulse_spec, spec_len,
                   tchns, base_chn, num_channels,
                   freq_num_fd, freqs_fd, y_ini);
       for(j=0;j<num_channels+2*freq_num_fd;j++)
       {
            double val = (y_ini[j]-y_mes[j])*(y_ini[j]-y_mes[j])*geo.w[j]*geo.w[j]/(num_channels+2*freq_num_fd);
            res += val;
       }
    }

    if(res>*residual) up[0]++;
    else up[0] = 0;
    *residual = res;
    for(i=0;i<nlay + 1;i++)
        x_ini[i] = x0[i];
        
    for(i=nlay+1;i<2*nlay;i++)
        dpth[i-nlay-1] = x0[i];
        
    for(i=0;i<POL_PAR_NUM*num_pol_lay;i++)
        cole_ini[i] = cole0[i];
        
    for(i = 0; i < num_channels+2*freq_num_fd+10; i++)
       free(Jacob[i]);

    free(Jacob);
    
    return 1;
}

