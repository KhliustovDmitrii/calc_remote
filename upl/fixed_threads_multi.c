/************************************\
*   AEM inversion with fixed net     *
*                                    *
*   Command line args:               *
*   file_in, file_out                *
*   [NTHREADS] default 1             *        
*   [pol_lay_num] if LAYERWISE       *
\************************************/

#include <stdio.h>
#include <complex.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <pthread.h>

#include "fft23.c"
#include "fdfun.c"
#include "math_utils.c"
#include "kalman.c"
#include "text_utils.c"
#include "uncert.c"

#include "globals.h"

/*----------------------------------------------*\

- Check for waveform file
- Check for taking fabs in frequency computations
- Check for ppm or nTl
- Check for polarization model in setting upper, lower and cole
- Check for fd channel order
- Check for data transform
\*----------------------------------------------*/

void get_configs(char *, char *, int);
void form_path_to_waveform(char *, char *, char *, int);
int parse(char *, double *);
int flinversion_fixed(geometry, int, double *, double *, int, int *, double *, double *, double *, int, double complex *,
                int, int *, int, int, int, double *, double *, int *, double *, double *, double *, double *, double *); 
void var_mes(double *, double *, double *, int);
void rho_mes(double, double *, double *, int);

int main(int argc, char **argv)
{
    printf("\n+++++++++++++++++++++++++ inversion started +++++++++++++++++++++++++\n");
    fflush(stdout);
    
    if(argc > 3)
      NTHREADS = atoi(argv[3]);
    else
      NTHREADS = 1;
    
    char buf[4096], tmp[256];

    int i, j, p;
    size_t n = 0;
    double *values, *args, *freqs;

    values = (double *)malloc(100*sizeof(double));
    memset(buf, 0, sizeof(buf));
    
    get_configs(argv[0], argv[1], PROG_REGIME);
    
    if(LAYERWISE == 1)
    {
      configs.pol_num = 1;
      configs.pol_inds[0] = atoi(argv[4]);
    }
   
    FILE *fwf;
    form_path_to_waveform(argv[0], buf, argv[1], PROG_REGIME);
    fwf = fopen(buf, "rt");
    double waveform[configs.spec_len];
    for(i=0;fgets(buf,100,fwf);i++)
        waveform[i] = atof(buf);


    printf("\n+++++++++++++++++++++++++ conf params obtained +++++++++++++++++++++++++\n");
    printf("freq_num_fd               %d\n", configs.freq_num_fd);
    printf("freqs_fd:\n");
    for(i=0; i<configs.freq_num_fd; i++)
       printf("%lf ", configs.freqs_fd[i]);
    printf("\nlay_num                   %d\n", configs.lay_num);
    printf("pol_num                   %d\n", configs.pol_num);
    printf("pol_inds:\n");
    for(i=0; i<configs.pol_num; i++)
       printf("%d ", configs.pol_inds[i]);
    printf("\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
    fflush(stdout);

    int tmp_flag;
    char *tmp_chr;
    char *data;
    char time[configs.time_size];
    double y_mes[configs.num_channels+2*configs.freq_num_fd];
    double y_mes_single_freq[2];
    int half_space_freq_num = 0;
    double cum_depth;
    double rho[configs.lay_num];
    double cole[POL_PAR_NUM*configs.pol_num];
    double rho_DA_ini[4];
    double x_ini[2*configs.lay_num + 1];
    double x_start[2*configs.lay_num + 1];
    double cole_ini[POL_PAR_NUM*configs.pol_num];
    double S_ini[(configs.lay_num + 1 + POL_PAR_NUM*configs.pol_num)*(configs.lay_num+1 + POL_PAR_NUM*configs.pol_num)];
    double S_pre[(configs.lay_num + 1 + POL_PAR_NUM*configs.pol_num)*(configs.lay_num+1 + POL_PAR_NUM*configs.pol_num)];
    double S0[(configs.lay_num+1 + POL_PAR_NUM*configs.pol_num)*(configs.lay_num+1 + POL_PAR_NUM*configs.pol_num)];
    double y_ini[configs.num_channels+2*configs.freq_num_fd];
    double upper[configs.lay_num + 1 + POL_PAR_NUM*configs.pol_num];
    double lower[configs.lay_num + 1 + POL_PAR_NUM*configs.pol_num];
    double upper_ini[2];
    double lower_ini[2];       
    double kde_lay_thick[configs.num_channels];
    double kde_res[configs.num_channels];
    int up = 0;
    int s7 = AVERAGE, s7c = 0;
    double mesv[configs.num_channels+2*configs.freq_num_fd];
    double res = -1;
    double v, v1;
    double rho_ini = RES_INI;
    int empty_pol_inds[1];
    double empty_cole[1];
    empty_pol_inds[0] = -1;
    empty_cole[0] = 1;
    double alta = 0;
    double vda = 0;
    int ft = 1;
    int start_with_rho_ini = 1;
    int data_cntr = 0;
    double dpth[configs.lay_num],d=configs.first_thick;
    int nlay = configs.lay_num;
    double weight = 1000.;
    int itr;
    double d_ini[4];
    double reg_time;
    double mes_buf[configs.lay_num + 1 + POL_PAR_NUM*configs.pol_num];
    double mes_buf_sup[configs.lay_num + 1 + POL_PAR_NUM*configs.pol_num];
    for(i=0;i<configs.lay_num;i++)
        upper[i] = MAX_RES;
    upper[configs.lay_num] = 5;
    for(i=0;i<configs.lay_num;i++)
        lower[i] = MIN_RES;
    lower[configs.lay_num] = -5;
    
    for(i = 0; i < configs.pol_num; i++)
    {
       switch(POL_MODEL)
        {
        case 0: //Cole-Cole   
           upper[POL_PAR_NUM*i + configs.lay_num + 1] = MAX_RES;
           upper[POL_PAR_NUM*i + configs.lay_num + 2] = 1;
           upper[POL_PAR_NUM*i + configs.lay_num + 3] = 1;
           lower[POL_PAR_NUM*i + configs.lay_num + 1] = MIN_RES;
           lower[POL_PAR_NUM*i + configs.lay_num + 2] = 0.000001;
           lower[POL_PAR_NUM*i + configs.lay_num + 3] = 0.000001;
           break;
        case 1: //Hyperbolic
           upper[POL_PAR_NUM*i + configs.lay_num + 1] = 1000;
           upper[POL_PAR_NUM*i + configs.lay_num + 2] = 1000;
           lower[POL_PAR_NUM*i + configs.lay_num + 1] = -1000;
           lower[POL_PAR_NUM*i + configs.lay_num + 2] = -1000;  
           break;
        }
    }
    
    double y_ini_full[configs.spec_len];
    double complex impulse_spec[configs.spec_len];
    
    FFT imp_fft;
    
    memset(&imp_fft, 0, sizeof(imp_fft));
    
    init_fft(&imp_fft,configs.spec_len);
    
    for(i=0;i<configs.spec_len;i++)
        imp_fft.xn[i] = waveform[i];
        
    fft_pro(&imp_fft,0);
    
    for(i=0;i<configs.spec_len;i++)
        impulse_spec[i] = imp_fft.fn[i];

    memset(mesv,0,sizeof(mesv));

    memset(time, 0, sizeof(time));

    d = configs.first_thick;
    for(i=0;i<configs.lay_num-1;i++,d*=configs.step) dpth[i] = d;

    data = buf + configs.time_size; 

    FILE *fin  = fopen(argv[1],"rt");
    FILE *fout = fopen(argv[2],"wt");
    
    printf("\n+++++++++++++++++++++++++  data files opened  +++++++++++++++++++++++++\n");
    fflush(stdout);

    double mom = 5050000000;

    geo.prim = 1./mom; 
    geo1.prim = 1./mom;
    geo_ini.prim = 1./mom;
        
    for(i = 0;i<configs.pol_num;i++)
    {
        switch(POL_MODEL)
        {
        case 0: //Cole-Cole   
           cole[POL_PAR_NUM*i] = RES_INI;
           cole[POL_PAR_NUM*i+1] = 0.001;
           cole[POL_PAR_NUM*i+2] = 0.5;
           break;
        case 1: //Hyperbolic
           cole[POL_PAR_NUM*i] = 0.1;
           cole[POL_PAR_NUM*i+1] = 1;   
           break;
        }
    }

    memset(x_ini,0,sizeof(x_ini));
    memset(cole_ini,0,sizeof(cole_ini));
    memset(S_ini,0,sizeof(S_ini));
    memset(buf,0,sizeof(buf));
    x_ini[configs.lay_num] = 0;
    for(i=0;i<POL_PAR_NUM*configs.pol_num;i++)
        cole_ini[i] = cole[i];

    for(i=0;fgets(buf,4096,fin);i++) 
    {
        if(buf[0]=='/') 
        {  
            memset(buf,0,sizeof(buf));
            continue;
        }

        if(buf[0]=='L' || buf[0]=='B' || buf[0]=='T') 
        { 
            fputs(buf,fout);
            printf("%s",buf);

            for(i = 0;i<configs.pol_num;i++)
            {
               switch(POL_MODEL)
               {
               case 0: //Cole-Cole   
                  cole[POL_PAR_NUM*i] = RES_INI;
                  cole[POL_PAR_NUM*i+1] = 0.001;
                  cole[POL_PAR_NUM*i+2] = 0.5;
                  break;
               case 1: //Hyperbolic
                  cole[POL_PAR_NUM*i] = 0.1;
                  cole[POL_PAR_NUM*i+1] = 1;   
                  break;
               }
            }
            memset(buf,0,sizeof(buf));
            ft = 1;
            continue;
        }

        memcpy(time,buf,sizeof(time)-1); // reading time

        if(strstr(data,"*")) 
        {
            memset(buf,0,sizeof(buf));
            continue;
        }
        parse(data, values);
        geo.hor_dist = values[configs.hor_dist_pos];
        geo.ver_dist = values[configs.ver_dist_pos];
        geo.alt = values[configs.alt_pos];

        for(i = 0; i<configs.freq_num_fd; i++)
        {
           if(FD_COMP_ORDER==0)
           {
              y_mes[2*i] = values[configs.first_mes_pos + configs.freq_num_fd + i];
              y_mes[2*i + 1] = values[configs.first_mes_pos + i]; 
           }
           else
           {
              y_mes[2*i] = values[configs.first_mes_pos + i];
              y_mes[2*i + 1] = values[configs.first_mes_pos + configs.freq_num_fd + i]; 
           } 
        }
        for(i = 0; i<configs.num_channels; i++)
           y_mes[i+2*configs.freq_num_fd] = values[configs.first_chan_pos + i];

        memset(buf,0,sizeof(buf));
        data_cntr++;

        for(int i=0;i<2*configs.freq_num_fd+configs.num_channels;i++)
            mesv[i] += y_mes[i];
            
        alta += geo.alt;
        vda += geo.ver_dist;
        s7c++;
        if(--s7)
            continue;
            
            
        for(int i=0;i<2*configs.freq_num_fd+configs.num_channels;i++)
            y_mes[i] = mesv[i]/s7c;
        geo.alt = alta/s7c;
        geo.ver_dist = vda/s7c;
        s7c = 0;
        s7 = AVERAGE;
        memset(mesv,0,sizeof(mesv));
        alta = 0;
        vda = 0;
        
        //Use difference of inphase component
        for(i=0;i<configs.freq_num_fd;i++) 
        {
           if(i<configs.freq_num_fd-1)
              y_mes[2*i+1] = y_mes[2*(i+1)+1] - y_mes[2*i+1];
           else
           {
             if(configs.freq_num_fd > 1)
               y_mes[2*i+1] = y_mes[2*(i-1)+1];
           }
        }

        //Take absolute values
        if(TAKE_ABS==1)
         for(i=0;i<2*configs.freq_num_fd;i++)
            y_mes[i] = fabs(y_mes[i]);
        else
         for(i=0;i<configs.freq_num_fd;i++)
            y_mes[2*i + 1] = -y_mes[2*i + 1];
        
        res = -1;
        up = 0;
        itr = 0;
        upper_ini[0] = MAX_RES;
        upper_ini[1] = 0;
        lower_ini[0] = MIN_RES;
        lower_ini[1] = 0;
        switch(INITIAL_METHOD)
        {
        case 0: //Start with manually set hs
           for(i = 0; i < configs.lay_num; i++) x_start[i] = RES_INI;
           rho_ini = RES_INI;
           break;
           
        case 1: //Unfiorm half-space computed for lowest frequency
           y_mes_single_freq[0] = y_mes[2*half_space_freq_num];
           y_mes_single_freq[1] = y_mes[2*half_space_freq_num + 1];
           rho_DA_ini[0] = RES_INI;
           rho_DA_ini[1] = 0;
           rho_DA_ini[2] = 0;
           rho_DA_ini[3] = 0;
           for (itr = 0;itr < MAX_ITER; itr++) 
           {
               memset(d_ini, 0, sizeof(d_ini));
               d_ini[0] = configs.ERR_INI;
               d_ini[3] = 0.1;
               up = 0;
               
               flinversion_fixed(geo,1,rho_DA_ini,empty_cole, 0, empty_pol_inds,dpth,y_ini,y_ini_full,configs.freq_num_td,impulse_spec,configs.spec_len,configs.tchns,configs.base_chn,0,1,y_mes_single_freq,
               &res,&up, d_ini, &(configs.freqs_fd[half_space_freq_num]), configs.freqs_td, upper_ini, lower_ini);
               
               rho_ini = rho_DA_ini[0];
               if(sqrt(res) <STOP_VAL) break;
               if(up) break;
           }
           for(i = 0; i < configs.lay_num; i++) x_start[i] = rho_ini;
           break;
           
        case 2: //Uniform half-space computed for lowest frequency with signal > 2 sigma
           geo_ini.ver_dist = geo.ver_dist;
           geo_ini.hor_dist = geo.hor_dist;
           geo_ini.alt = geo.alt;
           i = 0;
           while((y_mes[2*i] < 2/geo.w[2*i] || y_mes[2*i+1] < 2/geo.w[2*i+1]) && i < configs.freq_num_fd-2)
              i++;
              
           printf("\n For hs calculation use freq #%d\n", i);
              
           y_mes_single_freq[0] = y_mes[2*i];
           y_mes_single_freq[1] = y_mes[2*i + 1];
           geo_ini.w[0] = geo.w[2*i];
           geo_ini.w[1] = geo.w[2*i+1];
           rho_DA_ini[0] = RES_INI;
           rho_DA_ini[1] = 0;
           rho_DA_ini[2] = 0;
           rho_DA_ini[3] = 0;
           for (itr = 0;itr < MAX_ITER; itr++) 
           {
               memset(d_ini, 0, sizeof(d_ini));
               d_ini[0] = configs.ERR_INI;
               d_ini[3] = 0.1;
               up = 0;
               
               flinversion_fixed(geo_ini,1,rho_DA_ini,empty_cole, 0, empty_pol_inds,dpth,y_ini,y_ini_full,configs.freq_num_td,impulse_spec,configs.spec_len,configs.tchns,configs.base_chn,0,1,y_mes_single_freq,
               &res,&up, d_ini, &(configs.freqs_fd[i]), configs.freqs_td, upper_ini, lower_ini);
               
               rho_ini = rho_DA_ini[0];
               if(sqrt(res) <STOP_VAL) break;
               if(up) break;
           }
           for(i = 0; i < configs.lay_num; i++) x_start[i] = rho_ini;
           break;
        
        case 3: // KDE
           if(configs.num_channels==0) //No time channels to calculate
              for(i = 0; i < configs.lay_num; i++) x_start[i] = RES_INI;
           else //FIt HS and use Trigubovich' formula for skin layer thickness
           {
              geo_ini.ver_dist = geo.ver_dist;
              geo_ini.hor_dist = geo.hor_dist;
              geo_ini.alt = geo.alt;
              for(i = 0; i < configs.num_channels; i++)
              {
                  rho_DA_ini[0] = RES_INI;
                  rho_DA_ini[1] = 0;
                  rho_DA_ini[2] = 0;
                  rho_DA_ini[3] = 0;
                  geo_ini.w[0] = geo.w[2*configs.freq_num_fd + i];  
                  
                  for (itr = 0;itr < MAX_ITER; itr++)
                  {
                      memset(d_ini, 0, sizeof(d_ini));
                      d_ini[0] = configs.ERR_INI;
                      d_ini[3] = 0.1;
                      up = 0;
                     
                      flinversion_fixed(geo_ini, 1, rho_DA_ini, empty_cole,0,empty_pol_inds,dpth, &(y_ini[2*configs.freq_num_fd+i]), y_ini_full, configs.freq_num_td, impulse_spec,
                              configs.spec_len, &(configs.tchns[i]), configs.base_chn, 1, 0, &(y_mes[2*configs.freq_num_fd + i]),
                              &res,&up, d_ini, configs.freqs_fd, configs.freqs_td, upper_ini, lower_ini);

                      if(sqrt(res) <STOP_VAL) break;
                  }
                  
                  kde_res[i] = rho_DA_ini[0];
                  reg_time = 12.5*(configs.tchns[i] + configs.tchns[i+1])/2; //5 mks (1 in tchns) and mu0 are cancelled
                  kde_lay_thick[i] = sqrt(kde_res[i]*reg_time);
              }
              
              printf("KDE parameters \n");
              printf("\n-------------------------------------------------------------------------------------\n");
              printf("RESISTIVITIES:\n");
              for(i=0; i<configs.num_channels; i++)
                 printf("%f ",kde_res[i]);
              printf("\n");
              printf("THICKNESSES:\n");
              for(i=0; i<configs.num_channels; i++)
                 printf("%f ",kde_lay_thick[i]);
              printf("\n");
              
              //Discretize skin layers into fixed depth net
              cum_depth = dpth[0] + dpth[1];
              
              j = 0;
              x_start[0] = kde_res[0];
              for(i = 1; i<configs.lay_num; i++) //Fill all resistivities
              {
                  if(kde_lay_thick[j] > cum_depth && j < configs.num_channels - 1) // If above border
                  {
                     x_start[i] = kde_res[j]; // Set resistivity
                     cum_depth += dpth[i];
                  }
                  else // If crossed the border
                  {
                     if(j < configs.num_channels - 1)  // If not top KDE layer
                     {
                        j++; // Go to next KDE layer
                        i--; // And try to set current resistivity again
                     }
                     else // If top KDE layer
                     {
                        i--;
                        break; //Stop the cycle...
                     }
                  }
              }
              for(; i<configs.lay_num; i++) //... and fill remaining resisitivities with top KDE resistivity
                 x_start[i] = kde_res[configs.num_channels-1];
           }
           break;
        }



        res = sqrt(res);
       
        int iter = 0;
        x_ini[configs.lay_num] = 0;
        double Sr[(configs.lay_num + 1 + POL_PAR_NUM*configs.pol_num)*(configs.lay_num+1+POL_PAR_NUM*configs.pol_num)];
        res = -1;
        up = 0;

        memset(S_ini,0,sizeof(S_ini));
        for(i=configs.lay_num-1;i>=0;i--) 
        {
            if(i==nlay-1)
               S_ini[i+(configs.lay_num+1+POL_PAR_NUM*configs.pol_num)*i] = configs.ERR_INI;
            else 
            {
               S_ini[i+1+(configs.lay_num+1+POL_PAR_NUM*configs.pol_num)*i] = configs.COR_INI*configs.ERR_INI/S_ini[i+1+(configs.lay_num+1+POL_PAR_NUM*configs.pol_num)*(i+1)];
               S_ini[i+(configs.lay_num+1+POL_PAR_NUM*configs.pol_num)*i] = sqrt(configs.ERR_INI*configs.ERR_INI-S_ini[i+1+(configs.lay_num+1+POL_PAR_NUM*configs.pol_num)*i]*S_ini[i+1+(configs.lay_num+1+POL_PAR_NUM*configs.pol_num)*i]);
            }

        }
        S_ini[configs.lay_num+ configs.lay_num*(configs.lay_num+1+POL_PAR_NUM*configs.pol_num)] = 0.25;
        for(i = configs.lay_num+1; i < configs.lay_num + 1 + POL_PAR_NUM*configs.pol_num; i++) S_ini[i + i*(configs.lay_num+1+POL_PAR_NUM*configs.pol_num)] = 0.1;
        
        if(ft) 
        {
            memset(x_ini, 0, sizeof(x_ini));
            for(i=configs.lay_num-1;i>=0;i--)   
                x_ini[i] = x_start[i];
            ft = !ft;
            memcpy(S0,S_ini,sizeof(Sr));
        } 
        else 
        {
            for(i=configs.lay_num-1;i>=0;i--)   
                x_ini[i] = (x_start[i]*(1-1./weight)+x_ini[i]/weight);
                
            x_ini[configs.lay_num] = 0;
            
            for(i = 0; i < POL_PAR_NUM*configs.pol_num; i++) cole_ini[i] = cole[i]*(1-1./weight)+cole_ini[i]/weight;
            if(POL_MODEL==0)
               for(j = 0; j < configs.pol_num; j++) cole_ini[POL_PAR_NUM*j] = x_ini[configs.pol_inds[j]];
        }
        memcpy(Sr,S_ini,sizeof(Sr));
        
        //Print initial conditions
        printf("INITIAL CONDITIONS \n");
        printf("\n-------------------------------------------------------------------------------------\n");
        printf("TIME: %s \n ", time);
        printf("NUM LAYERS:    %d\n", configs.lay_num);
        printf("RESISTIVITIES:\n");
        for(i=0; i<configs.lay_num; i++)
           printf("%f ",x_ini[i]);
        if(configs.pol_num > 0)
        {
           printf("\nCOLE:\n");
           for(i=0; i<POL_PAR_NUM*configs.pol_num; i++)
              printf("%f ",cole_ini[i]);
        }
        printf("\n-------------------------------------------------------------------------------------\n");      

        
        //Produce inversion

        memset(mes_buf_sup, 0, sizeof(mes_buf_sup));
        for (iter = 0;iter < MAX_ITER; iter++) 
        {
            memcpy(S_ini,Sr,sizeof(Sr));
            
            if(TD_ONLY)
            {
               geo1.ver_dist = geo.ver_dist;
               geo1.hor_dist = geo.hor_dist;
               geo1.alt = geo.alt;
               flinversion_fixed(geo1, configs.lay_num, x_ini, cole_ini,configs.pol_num,configs.pol_inds,dpth, &(y_ini[2*configs.freq_num_fd]), y_ini_full, configs.freq_num_td, impulse_spec,
               configs.spec_len, configs.tchns, configs.base_chn, configs.num_channels, 0, &(y_mes[2*configs.freq_num_fd]),
               &res,&up, S_ini, configs.freqs_fd, configs.freqs_td, upper, lower);
            }
            else
               flinversion_fixed(geo, configs.lay_num, x_ini, cole_ini,configs.pol_num,configs.pol_inds,dpth, y_ini, y_ini_full, configs.freq_num_td, impulse_spec,
               configs.spec_len, configs.tchns, configs.base_chn, configs.num_channels, configs.freq_num_fd, y_mes,
               &res,&up, S_ini, configs.freqs_fd, configs.freqs_td, upper, lower);
            
            var_mes(S_ini, Sr, mes_buf, configs.lay_num + 1 + POL_PAR_NUM*configs.pol_num);
            for(i = 0; i < configs.lay_num + 1 + POL_PAR_NUM*configs.pol_num; i++)
               if(mes_buf[i] > mes_buf_sup[i])
                  mes_buf_sup[i] = mes_buf[i];
            
            
            if(sqrt(res) <STOP_VAL) break;
            if(up) break;
        }
        
        res = sqrt(res);
        weight = (res>1)?res:1.;

        fprintf(fout, "%s %f %d     ",time, res, iter);
        
        p = 0;   
        switch(POL_MODEL)
        {
        case 0:
           for(i=0; i<configs.lay_num+1; i++)
              fprintf(fout, "%f ",x_ini[i]);
           break;
        case 1:
           for(i=0; i<configs.lay_num+1; i++)
           {
              if(i<configs.lay_num && configs.pol_inds[p] == i)
              {
                fprintf(fout, "%f ",x_ini[i] + 1/cole_ini[POL_PAR_NUM*p + 1]);
                p++;
              }
              else
                fprintf(fout, "%f ",x_ini[i]);                
           }
           break;
        }
        
        double depp = 0;
        for(i=0;i<configs.lay_num-1;i++) 
        {
            fprintf(fout,"%.3f ",depp+dpth[i]*.5);
            depp += dpth[i];
        }

        fprintf(fout,"%.3f %.3f ",depp+dpth[configs.lay_num-2]*.5, rho_ini);
        
        for(i=0;i<POL_PAR_NUM*configs.pol_num;i++) 
           fprintf(fout,"%f ",cole_ini[i]);
           
        for(i = 0; i < configs.lay_num + 1 + POL_PAR_NUM*configs.pol_num; i++) //Print sup estimability measures
           fprintf(fout, "%f ", mes_buf_sup[i]);
           
           
        if(VERBOSE>0) //Print final step estimability measures
           for(i = 0; i < configs.lay_num + 1 + POL_PAR_NUM*configs.pol_num; i++)
              fprintf(fout, "%f ", mes_buf[i]);
        
        if(INITIAL_METHOD<3 && VERBOSE>1) //Print rho estimability measures
        {
           rho_mes(rho_ini, x_ini, mes_buf, configs.lay_num);
           for(i = 0; i < configs.lay_num; i++)
              fprintf(fout, "%f ", mes_buf[i]);
        }

        fprintf(fout,"\n");
        fflush(fout);
        
        //Print output
        printf("---------------> %f <-----  %d \n", res, iter);
        printf("\n");
        printf("RESISTIVITIES:\n");
        for(i=0; i<configs.lay_num; i++)
           printf("%f ",x_ini[i]);
        printf("\n\nAltitude correction:\n");
        printf("----> %f <----\n",x_ini[configs.lay_num]);
                
        if(configs.pol_num>0)
        {
           printf("\nCOLE:\n");
           for(i=0; i<POL_PAR_NUM*configs.pol_num; i++)
              printf("%f ",cole_ini[i]);
        }
        printf("\n");
        printf("\n##################################################################################################\n");   
        
        fflush(stdout);
        
        continue;
    }
    printf("\n end");
    return 0;
}


