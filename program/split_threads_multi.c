/************************************\
*   AEM inversion with fixed net     *
*                                    *
*   Command line args:               *
*   file_in, file_out                *
*   [NTHREADS] default 1             *
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
int flinversion_free(geometry, int, double *, double *, int, int *, double *, double *, double *, int, double complex *,
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
    
    char buf[BUF_SIZE], tmp[256];

    int i, j, read_num;
    size_t n = 0;
    double *values, *args, *freqs;

    values = (double *)malloc(100*sizeof(double));

    memset(buf, 0, sizeof(buf));
    
    get_configs(argv[0], argv[1], PROG_REGIME);
     
    FILE *fwf;
    form_path_to_waveform(argv[0], buf, argv[1], PROG_REGIME);
    fwf = fopen(buf, "rt");
    double waveform[configs.spec_len];
    for(i=0;fgets(buf,100,fwf);i++)
        waveform[i] = atof(buf);

    printf("\n+++++++++++++++++++++++++ conf params obtained +++++++++++++++++++++++++\n");
    printf("freq_num_fd               %d\n", configs.freq_num_fd);
    printf("pol_inds:\n");
    for(i=0; i<configs.freq_num_fd; i++)
       printf("%lf ", configs.freqs_fd[i]);
    printf("pol_num                   %d\n", configs.pol_num);
    printf("pol_inds:\n");
    for(i=0; i<configs.pol_num; i++)
       printf("%d ", configs.pol_inds[i]);
    printf("\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
    fflush(stdout);

    char data[BUF_SIZE];
    char time_val[configs.time_size+5];
    double y_mes[configs.num_channels+2*configs.freq_num_fd];
    double y_mes_single_freq[2];
    int half_space_freq_num = 0;
    double rho[configs.lay_num];
    double cole[POL_PAR_NUM*configs.pol_num];
    double cole_best[POL_PAR_NUM*configs.pol_num];
    double rho_DA_ini[4];
    double x_ini[2*MAX_FREE_LAYERS], x_temp[2*MAX_FREE_LAYERS];
    double x_best_loc[2*MAX_FREE_LAYERS];
    double cole_best_loc[POL_PAR_NUM*configs.pol_num];
    double x_best[2*MAX_FREE_LAYERS];
    double cole_ini[POL_PAR_NUM*configs.pol_num];
    double S_ini[(2*MAX_FREE_LAYERS + POL_PAR_NUM*configs.pol_num)*(2*MAX_FREE_LAYERS + POL_PAR_NUM*configs.pol_num)];
    double S_pre[(2*MAX_FREE_LAYERS + POL_PAR_NUM*configs.pol_num)*(2*MAX_FREE_LAYERS + POL_PAR_NUM*configs.pol_num)];
    double S0[(2*MAX_FREE_LAYERS + POL_PAR_NUM*configs.pol_num)*(2*MAX_FREE_LAYERS + POL_PAR_NUM*configs.pol_num)];
    double S[(2*MAX_FREE_LAYERS + POL_PAR_NUM*configs.pol_num)*(2*MAX_FREE_LAYERS + POL_PAR_NUM*configs.pol_num)];
    double Sr[(2*MAX_FREE_LAYERS + POL_PAR_NUM*configs.pol_num)*(2*MAX_FREE_LAYERS + POL_PAR_NUM*configs.pol_num)];
    double Sr_best[(2*MAX_FREE_LAYERS + POL_PAR_NUM*configs.pol_num)*(2*MAX_FREE_LAYERS + POL_PAR_NUM*configs.pol_num)];
    double Sr_best_loc[(2*MAX_FREE_LAYERS + POL_PAR_NUM*configs.pol_num)*(2*MAX_FREE_LAYERS + POL_PAR_NUM*configs.pol_num)];
    double Sr_temp[(2*MAX_FREE_LAYERS + POL_PAR_NUM*configs.pol_num)*(2*MAX_FREE_LAYERS + POL_PAR_NUM*configs.pol_num)];
    double S_best[(2*MAX_FREE_LAYERS + POL_PAR_NUM*configs.pol_num)*(2*MAX_FREE_LAYERS + POL_PAR_NUM*configs.pol_num)];
    double y_ini[configs.num_channels+2*configs.freq_num_fd];
    double upper[2*MAX_FREE_LAYERS + POL_PAR_NUM*configs.pol_num];
    double lower[2*MAX_FREE_LAYERS + POL_PAR_NUM*configs.pol_num];
    double upper_ini[2];
    double lower_ini[2];        
    int up = 0;
    int s7 = AVERAGE, s7c = 0;
    double mesv[configs.num_channels+2*configs.freq_num_fd];
    double alta = 0;
    double vda = 0;
    int ft = 1;
    int start_with_rho_ini = 1;
    int data_cntr = 0;
    double dpth[MAX_FREE_LAYERS], dpth_temp[MAX_FREE_LAYERS], d=configs.first_thick;
    double dpth_best[MAX_FREE_LAYERS];
    double dpth_best_loc[MAX_FREE_LAYERS];
    double weight = 1000.;
    double v, v1;
    double mes_buf[2*configs.lay_num + 1 + POL_PAR_NUM*configs.pol_num];
    
    
    
    
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

    FILE *fin  = fopen(argv[1],"rt");
    FILE *fout = fopen(argv[2],"wt");
    printf("data files opened\n");
    fflush(stdout);
    
    double mom = 5050000000;
    geo.prim = 1./mom;                 
    geo1.prim = 1./mom;

    memset(x_ini,0,sizeof(x_ini));
    memset(cole_ini,0,sizeof(cole_ini));
    memset(S_ini,0,sizeof(S_ini));

    for(i=0;fgets(buf,BUF_SIZE,fin);i++) 
    {
        if(buf[0]=='/')  // reading comments
            continue;

        if(buf[0]=='L' || buf[0]=='B' || buf[0]=='T') // reading Lines or Base or Trends
        {
            fputs(buf,fout);
            printf("%s",buf);
            ft = 1;
            continue;
        }

        strncpy(data, &(buf[configs.time_size]), BUF_SIZE-128-configs.time_size);

        if(strstr(data,"*")) // all the nonmeasurements are skipped
            continue;
        
        parse(data, values);
        ft = 1; 
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

        data_cntr++;

    // averaging for AVERAGE samples
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
        
        geo1.ver_dist = geo.ver_dist;
        geo1.hor_dist = geo.hor_dist;
        geo1.alt = geo.alt;

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
        

//########################################################################################################
    
    int num_free_layers;
    int num_in_gen = 0;
    int gen_pow;
    int p;
    int iter;
    int lays_to_split;
    int last_update = 1; //Counter for layer where we updated solution
    double res1 = 10000000;
    double res_loc;
    double res = -1;
    int stop_flag = 0;
    int final_lay_num;
    int fin_pol_num;
    memset(dpth, 0, sizeof(dpth));
    memset(dpth_temp, 0, sizeof(dpth_temp));
    memset(x_ini, 0, sizeof(x_ini));
    memset(y_ini, 0, sizeof(y_ini));
    memset(x_temp, 0, sizeof(x_temp));
    memset(x_best, 0, sizeof(x_best));
    memset(dpth_best, 0, sizeof(dpth_best));
    memset(x_best_loc, 0, sizeof(x_best_loc));
    memset(dpth_best_loc, 0, sizeof(dpth_best_loc));
    memset(S_best, 0, sizeof(S_best));
    memset(S, 0, sizeof(S));
    memset(Sr_best, 0, sizeof(Sr_best));
    memset(Sr_best_loc, 0, sizeof(Sr_best_loc));
    memset(Sr, 0, sizeof(Sr));
    memset(Sr_temp, 0, sizeof(Sr_temp));
    memset(cole_ini, 0, sizeof(cole_ini));
    memset(cole_best, 0, sizeof(cole_best));
    memset(cole_best_loc, 0, sizeof(cole_best_loc));
    int loc_pol_inds[configs.pol_num];
    int loc_pol_num;
    
    res1 = 10000000;
    for(num_free_layers = 1; num_free_layers < MAX_FREE_LAYERS+1; num_free_layers++)
    {
       
       if(stop_flag==0)
       {
        
           if(num_free_layers == 1) //Define number of layers to be splitted
              lays_to_split = 1;
           else
              lays_to_split = num_free_layers-1;
           

           memcpy(x_temp, x_best_loc, sizeof(x_best_loc)); //Use previous best solution as starting model
           memcpy(dpth_temp, dpth_best_loc, sizeof(dpth_best_loc));
           memcpy(Sr_temp, Sr, sizeof(Sr));
           
           for(num_in_gen = 0; num_in_gen < lays_to_split; num_in_gen++)//Iterate over all layers to split
           {
              //Form initial conditions
              if(num_free_layers==1) //If half-space
              {
                 x_ini[0] = RES_INI; 
                 x_ini[1] = 0;
              }
              else
              {
                  //Copy resistivities and depths for layers we do not split
                  memcpy(x_ini, x_temp, sizeof(x_temp));
                  memcpy(dpth, dpth_temp, sizeof(dpth_temp));
                  
                  //Split the layer
                  if(num_in_gen==num_free_layers-2) //If the last layer in generation
                  {
                     x_ini[num_in_gen+1] = x_temp[num_in_gen];
                     dpth[num_in_gen] = 30;
                  }
                  else //Otherwise
                  {
                     x_ini[num_in_gen+1] = x_temp[num_in_gen];
                     dpth[num_in_gen+1] = 0.5*dpth_temp[num_in_gen];
                     dpth[num_in_gen] = 0.5*dpth_temp[num_in_gen];
                  }
                  //Copy resistivities and depths for layers we do not split
                  for(i = num_in_gen+2; i < num_free_layers; i++)
                  {
                      x_ini[i] = x_temp[i-1];
                      dpth[i] = dpth_temp[i-1];
                  }
                  
                  //Zero altitude correction
                  x_ini[num_free_layers] = 0;   
              }   
              
              loc_pol_num = 0;
              memset(loc_pol_inds, 0, sizeof(loc_pol_inds));
              
              
              //Fill in polarization parameters 
              if(configs.pol_num>0)
              {
                 for(i = 0; i < configs.pol_num; i++)
                 {
                    if(configs.pol_inds[i]<=num_free_layers-1)
                    {
                       loc_pol_inds[i] = configs.pol_inds[i];
                       loc_pol_num++;
                       
                       switch(POL_MODEL)
                       {
                       case 0: //Cole-Cole   
                          cole_ini[POL_PAR_NUM*i] = 1.1*x_ini[configs.pol_inds[i]];
                          cole_ini[POL_PAR_NUM*i+1] = 0.001;
                          cole_ini[POL_PAR_NUM*i+2] = 0.5;
                          break;
                       case 1: //Hyperbolic
                          cole_ini[POL_PAR_NUM*i] = 0.1;
                          cole_ini[POL_PAR_NUM*i+1] = 1;   
                          break;
                       }
                    }
                 } 
              }

              //Fill in dispersions
              memset(S, 0, sizeof(S));
              memset(Sr, 0, sizeof(Sr));  
              
              //Print initial conditions
                 printf("INITIAL CONDITIONS \n");
                 printf("------------------------------------------------------------------------\n");
                 printf("NUM LAYERS:    %d\n", num_free_layers);
                 printf("\n");
                    
                 printf("RESISTIVITIES:\n");
                 switch(POL_MODEL)
                 {
                 case 0:
                    for(i=0; i<num_free_layers; i++)
                       printf("%f ",x_ini[i]);
                    break;
                 case 1:
                    p = 0;
                    for(i=0; i<num_free_layers; i++)
                    {
                       if(i<loc_pol_num && loc_pol_inds[p] == i)
                       {
                          printf("%f ",x_ini[i] + 1/cole_ini[POL_PAR_NUM*p + 1]);
                          p++;
                       }
                       else
                          printf("%f ",x_ini[i]);                
                    }
                   break;
                  }
                     
                  printf("\n\nAltitude correction:\n");
                    
                  printf("----> %f <----\n",x_ini[num_free_layers]);
                    
                  printf("\nDEPTHS:\n");
                  for(i=0; i<num_free_layers-1; i++)
                     printf("%f ",dpth[i]);
                    
                  printf("\nPOLARIZATION PARAMETERS:\n");
                  for(i=0; i<POL_PAR_NUM*loc_pol_num; i++)
                   printf("%f ",cole_ini[i]);
               
                 printf("\n------------------------------------------------------------------------\n");         
                 
              //Fill dispersions for resistivities
              for(i=num_free_layers-1; i>=0; i--) 
              {
                  if(i==num_free_layers-1)
                     S[i+(2*num_free_layers+POL_PAR_NUM*loc_pol_num)*i] = configs.ERR_INI;
                  else
                  {
                     S[i+1+(2*num_free_layers+POL_PAR_NUM*loc_pol_num)*i] = configs.COR_INI*configs.ERR_INI/S[i+1+(2*num_free_layers+POL_PAR_NUM*loc_pol_num)*(i+1)];
                     S[i+(2*num_free_layers+POL_PAR_NUM*loc_pol_num)*i] = sqrt(configs.ERR_INI*configs.ERR_INI-S[i+1+(2*num_free_layers+POL_PAR_NUM*loc_pol_num)*i]*S[i+1+(2*num_free_layers+POL_PAR_NUM*loc_pol_num)*i]);
                  }
              }
                 
             //Dispersion for altitude correction
             S[num_free_layers + num_free_layers*(2*num_free_layers+POL_PAR_NUM*loc_pol_num)] = 0.25;
                 
             //Dispersions for depths
             for(i=num_free_layers+1; i<2*num_free_layers; i++)
                 S[i + i*(2*num_free_layers+POL_PAR_NUM*loc_pol_num)] = 0.25;
                    
             //Cole-Cole dispersions
             for(i=2*num_free_layers; i<2*num_free_layers+POL_PAR_NUM*loc_pol_num; i++)
                 S[i + i*(2*num_free_layers+POL_PAR_NUM*loc_pol_num)] = 0.1;   
                 
             /*   
             printf("\nS:\n");
             for(i=0; i<2*num_free_layers + POL_PAR_NUM*loc_pol_num; i++)
             {
                for(j=0; j<2*num_free_layers + POL_PAR_NUM*loc_pol_num; j++)
                   printf("%f ",S[j + i*(2*num_free_layers + POL_PAR_NUM*loc_pol_num)]);
                printf("\n");
             }
             printf("\n##################################################\n");     
             */
             
             //Save dispersions for estimability measures calculation
                 
             memcpy(Sr, S, sizeof(Sr));
                 
             //Fll in parameter bounds
             for(i = 0; i < 2*num_free_layers; i++)
             {
                 if(i < num_free_layers)
                 {
                    lower[i] = MIN_RES;
                    upper[i] = MAX_RES;   
                 }
                 if(i == num_free_layers)
                 {
                    lower[i] = -5;
                    upper[i] = 5;   
                 }
                 if(i > num_free_layers && i < 2*num_free_layers)
                 {
                    lower[i] = 1;
                    upper[i] = 200;   
                 }
             }
             for(i = 0; i < loc_pol_num; i++)
             {
                switch(POL_MODEL)
                {
                case 0: //Cole-Cole   
                   upper[POL_PAR_NUM*i + 2*num_free_layers] = MAX_RES;
                   upper[POL_PAR_NUM*i + 1 + 2*num_free_layers] = 1;
                   upper[POL_PAR_NUM*i + 2 + 2*num_free_layers] = 1;
                   lower[POL_PAR_NUM*i + 2*num_free_layers] = MIN_RES;
                   lower[POL_PAR_NUM*i + 1 + 2*num_free_layers] = 0.000001;
                   lower[POL_PAR_NUM*i + 2 + 2*num_free_layers] = 0.000001;
                   break;
                case 1: //Hyperbolic
                   upper[POL_PAR_NUM*i + 2*num_free_layers] = 1000;
                   upper[POL_PAR_NUM*i + 1 + 2*num_free_layers] = 1000;
                   lower[POL_PAR_NUM*i + 2*num_free_layers] = -1000;
                   lower[POL_PAR_NUM*i + 1 + 2*num_free_layers] = -1000;   
                   break;
                }
             } 
              
             //Perform inversion
             for (iter = 0;iter < MAX_ITER; iter++)
             {
                 memcpy(S,Sr,sizeof(S));
                 
                 if(TD_ONLY)
                    flinversion_free(geo1, num_free_layers, x_ini, cole_ini, loc_pol_num, loc_pol_inds, dpth,
                                    &(y_ini[2*configs.freq_num_fd]), y_ini_full, configs.freq_num_td, impulse_spec, configs.spec_len, configs.tchns, configs.base_chn,
                                    configs.num_channels, 0, &(y_mes[2*configs.freq_num_fd]), &res, &up, S, configs.freqs_fd, configs.freqs_td,
                                    upper, lower);
                 else
                    flinversion_free(geo, num_free_layers, x_ini, cole_ini, loc_pol_num, loc_pol_inds, dpth,
                                    y_ini, y_ini_full, configs.freq_num_td, impulse_spec, configs.spec_len, configs.tchns, configs.base_chn,
                                    configs.num_channels, configs.freq_num_fd, y_mes, &res, &up, S, configs.freqs_fd, configs.freqs_td,
                                    upper, lower);
                                    
                 if(sqrt(res) <STOP_VAL) break;
             }
             res = sqrt(res);
                 
             //Print output
             printf("RESULT\n");
             printf("------------------------------------------------------------------------\n");
             printf("NUM LAYERS: %d  ---------------> %f <-----  %d \n", num_free_layers, res, iter);
             printf("\n");
                 
             printf("RESISTIVITIES:\n");
             switch(POL_MODEL)
             {
             case 0:
                for(i=0; i<num_free_layers; i++)
                    printf("%f ",x_ini[i]);
                break;
             case 1:
                p = 0;
                for(i=0; i<num_free_layers; i++)
                {
                    if(i<loc_pol_num && loc_pol_inds[p] == i)
                    {
                       printf("%f ",x_ini[i] + 1/cole_ini[POL_PAR_NUM*p + 1]);
                       p++;
                    }
                    else
                       printf("%f ",x_ini[i]);                
                }
                break;
             }
               
             printf("\n\nAltitude correction:\n");
             printf("----> %f <----\n",x_ini[num_free_layers]);
                 
             printf("\nDEPTHS:\n");
             for(i=0; i<num_free_layers-1; i++)
                printf("%f ",dpth[i]);
                
             printf("\nPOLARIZATION PARAMETERS:\n");
             for(i=0; i<POL_PAR_NUM*loc_pol_num; i++)
             {
                printf("%f ",cole_ini[i]);
             }
             
             printf("\n");
             printf("\n#########################################################################################################\n");   
                  
             if(num_in_gen==0 || res < res_loc)//If first model in generation or best model in generation
             {
                res_loc = res;
                memcpy(x_best_loc, x_ini, sizeof(x_ini));
                memcpy(cole_best_loc, cole_ini, sizeof(cole_ini));
                memcpy(dpth_best_loc, dpth, sizeof(dpth));   
                fin_pol_num = loc_pol_num;
             }
                     
             //If solution is better than previous best, save it
             if(res < res1)
             {
                 res1 = res;
                 memcpy(x_best, x_ini, sizeof(x_ini));
                 memcpy(dpth_best, dpth, sizeof(dpth));
                 memcpy(cole_best, cole_ini, sizeof(cole_ini));
                 memcpy(S_best, S, sizeof(S));
                 memcpy(Sr_best, Sr, sizeof(Sr));
                 last_update = num_free_layers;
                 final_lay_num = lays_to_split+1;
                 fin_pol_num = loc_pol_num;
             }
             //If solution is good enough or no updates for current number of layers
             if(res < STOP_VAL || ((num_in_gen == lays_to_split - 1 && last_update < num_free_layers) && (GREEDY>0)))
             {
                 printf("\n-------------->STOPPING INVERSION<-------------- \n");
                 stop_flag = 1; //Save it and stop calculating more layers
                 if(res < STOP_VAL)
                 {
                    memcpy(x_best, x_ini, sizeof(x_ini));
                    memcpy(dpth_best, dpth, sizeof(dpth));
                    memcpy(cole_best, cole_ini, sizeof(cole_ini));
                    memcpy(S_best, S, sizeof(S));
                    memcpy(Sr_best, Sr, sizeof(Sr));
                    final_lay_num = lays_to_split+1;
                    fin_pol_num = loc_pol_num;
                 }
                 else
                 {
                    memcpy(x_best, x_temp, sizeof(x_temp));
                    memcpy(dpth_best, dpth_temp, sizeof(dpth_temp));
                    memcpy(cole_best, cole_ini, sizeof(cole_ini));
                    memcpy(Sr_best, Sr_temp, sizeof(Sr_temp));
                    final_lay_num = lays_to_split;
                    fin_pol_num = loc_pol_num;
                 }
                 break;
             }
           }
        }
     }
     
     for(i = 0; i < configs.time_size; i++)
        fprintf(fout, "%c", buf[i]);
        
     fprintf(fout, " %f  %d   ", res1, final_lay_num);
        
     switch(POL_MODEL)
     {
     case 0:
        for(i=0; i<final_lay_num; i++)
           fprintf(fout, "%f ",x_best[i]);
        break;
     case 1:
        p = 0;
        for(i=0; i<final_lay_num; i++)
        {
           if(i<fin_pol_num && loc_pol_inds[p] == i)
           {
             fprintf(fout, "%f ",x_best[i] + 1/cole_best[POL_PAR_NUM*p + 1]);
             p++;
           }
           else
             fprintf(fout, "%f ",x_best[i]);                
        }
        break;
     }
     
     for(i=final_lay_num;i<MAX_FREE_LAYERS;i++) 
         fprintf(fout,"%.3f ",x_best[final_lay_num-1]);
     fprintf(fout,"%.3f ",x_best[final_lay_num]);

     double depp = 0;
     for(i=final_lay_num-1;i<MAX_FREE_LAYERS-1;i++)
        dpth_best[i] = 30;
        
     for(i=0;i<MAX_FREE_LAYERS-1;i++) 
     {
        fprintf(fout,"%.3f ",depp+dpth_best[i]*.5);
        depp += dpth_best[i];
     }
     fprintf(fout,"%.3f ",depp+dpth_best[MAX_FREE_LAYERS-2]*.5);
     
     var_mes(S_best, Sr_best, mes_buf, 2*final_lay_num+POL_PAR_NUM*fin_pol_num);
     
     for(i = 0; i < 2*final_lay_num+POL_PAR_NUM*fin_pol_num; i++)
         fprintf(fout, "%f ", mes_buf[i]);

         
      for(i=0; i<POL_PAR_NUM*fin_pol_num; i++)
         fprintf(fout, "%f ",cole_best[i]);

        fprintf(fout,"\n");
        fflush(fout);
    
    
    }
  
    printf("\n end");
    return 0;
}
