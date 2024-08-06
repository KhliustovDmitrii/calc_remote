/********************************************\
 * Routines for working with text data      *
\* ******************************************/

#include "globals.h"

int parse(char *buf, double *values)
{
    int is_digit = 0;
    int j = 0;
    int i;
    int read_num = 0;
    char tmp[256];
    
    for(i = 0; i < BUF_SIZE - 128; i++)
    {
        if((buf[i]>='0' && buf[i]<='9')||buf[i]=='.'||buf[i]=='-')
        {
           tmp[j] = buf[i];
           j++;
           is_digit = 1;
        } 
        else
        {
          if(is_digit==1)
          {
             values[read_num] = atof(tmp);
             read_num++;
             is_digit = 0;
             for(j = 0; j < 256; j++)
                tmp[j] = (char) 0;
          }
          j = 0;
        }
    }
    return 0;
}

void form_path_to_configs(char *prog, char *path, char *remote_path, int regime)
{
   char full_name[1024];
   strcpy(full_name, prog); 
   char *cur_name = strrchr(full_name, '/');
   char ququ[1024] = "configs"; 
   char path_to[2048];
   char *ptr_path_to;
   int slash_count = 0;
   int slash_i = 0;
   strcat(ququ, cur_name);
   strcat(ququ, ".conf");
 
   if(regime==0) //Local program
   {
      sprintf(path, "%s", ququ);
   }
   else //Server program
   {
      strcpy(path_to, remote_path);
 
      ptr_path_to = path_to;

      while(slash_count < 6 && path_to[slash_i]!='\0')
      {
        if(path_to[slash_i]=='/') slash_count++;
        slash_i++;
      }   
      path_to[slash_i] = '\0';

      strcat(path_to, ququ);
      
      sprintf(path, "%s", path_to);
   }
}

void form_path_to_waveform(char *prog, char *path, char *remote_path, int regime)
{
   char ququ[1024];
   char path_to[2048];
   char *ptr_path_to;
   int slash_count = 0;
   int slash_i = 0;
   if(regime==0) //Local program
   {
      sprintf(path, "%s", "configs/waveform.XYZ");
   }
   else //Server program
   {
      sprintf(ququ, "%s", "configs/waveform.XYZ"); 
      strcpy(path_to, remote_path);
 
      ptr_path_to = path_to;

      slash_count = 0;
      slash_i = 0;

      while(slash_count < 6 && path_to[slash_i]!='\0')
      {
         if(path_to[slash_i]=='/') slash_count++;
         slash_i++;
      }   
      path_to[slash_i] = '\0';

      strcat(path_to, ququ);
      sprintf(path, "%s", path_to);
   }
}

void get_configs(char *prog, char *remote_path, int regime)
{
   char path[1024], buf[BUF_SIZE];
   double values[100];
   int i;
   form_path_to_configs(prog, path, remote_path, regime);
   FILE *conf = fopen(path, "r");
   
   //Read number of fd freqs
   fgets(buf, 2000, conf);
   fgets(buf, 2000, conf);
   parse(buf, values);
   configs.freq_num_fd = (int)(values[0]+0.1);
    
//Read fd freqs
   fgets(buf, 2000, conf);
   fgets(buf, 2000, conf);
   parse(buf, values);
   for(i = 0; i < configs.freq_num_fd; i++)
      configs.freqs_fd[i] = values[i];

//Read lay_num, first_thick, step
   fgets(buf, 2000, conf);
   fgets(buf, 2000, conf);
   parse(buf, values);
    
   configs.lay_num = values[0];
   configs.first_thick = values[1];
   configs.step = values[2];

//Read ERR_INI, COR_INI
   fgets(buf, 2000, conf);
   fgets(buf, 2000, conf);
   parse(buf, values);

   configs.ERR_INI = values[0];
   configs.COR_INI = values[1];
   
//Read position of measurements, channels, hor_dist, ver_dist, alt
   fgets(buf, 2000, conf);
   fgets(buf, 2000, conf);
   parse(buf, values);

   configs.first_mes_pos = values[0];
   configs.last_mes_pos = values[1];
   configs.first_chan_pos = values[2];
   configs.last_chan_pos = values[3];
   configs.hor_dist_pos = values[4];
   configs.ver_dist_pos = values[5];
   configs.alt_pos = values[6];

//Read size of time string  
   fgets(buf, 2000, conf);
   fgets(buf, 2000, conf);
   parse(buf, values);

   configs.time_size = (int)(values[0]+0.1);
    
//Read number of pol layers
   fgets(buf, 2000, conf);
   fgets(buf, 2000, conf);
   parse(buf, values);
   
   configs.pol_num = (int)(values[0]+0.1);

//Read indices of pol layers  

   fgets(buf, 2000, conf);
   fgets(buf, 2000, conf);
   parse(buf, values);

   for(i = 0; i < configs.pol_num; i++)
    configs.pol_inds[i] = (int)(values[i]+0.1);
        
//Read base frequency, number of freqs, len of waveform, number of channels and base channel
   fgets(buf, 2000, conf);
   fgets(buf, 2000, conf);
   parse(buf, values);

   configs.base_freq = values[0];
   configs.freq_num_td = (int)(values[1]+0.1);
   configs.spec_len = (int)(values[2]+0.1);
   configs.base_chn = (int)(values[3]+0.1);
   configs.num_channels = (int)(values[4]+0.1);
    
//Read borders of time channels
    
   fgets(buf, 2000, conf);
   fgets(buf, 2000, conf);
   parse(buf, values);
    
   for(i = 0; i < configs.num_channels + 1; i++)
    configs.tchns[i] = (int)(values[i]+0.1);
    
   memset(configs.freqs_td, 0, sizeof(configs.freqs_td));
   for(i=1;i<2*configs.freq_num_td;i+=2)
      configs.freqs_td[(i-1)/2] = configs.base_freq*i;
    
    
//Read R_i_i
   fgets(buf, 2000, conf);
   fgets(buf, 2000, conf);
   parse(buf, values);
    
   geo.w = (double *)malloc((configs.num_channels+2*configs.freq_num_fd)*sizeof(double));
   geo1.w = (double *)malloc((configs.num_channels)*sizeof(double));
   geo_ini.w = (double *)malloc((configs.num_channels)*sizeof(double));
   
   memset(geo.w,0,sizeof(geo.w));
   memset(geo1.w,0,sizeof(geo1.w));
   memset(geo_ini.w,0,sizeof(geo1.w));
   for(i = 0; i < configs.freq_num_fd; i++)
   {
      geo.w[2*i] = 1/values[i];
      geo.w[2*i+1] = 1/values[i];
   }
    
   for(i=0;i<configs.freq_num_fd-1;i++)
      geo.w[i*2+1] = 1./sqrt(1./(geo.w[i*2+1]*geo.w[i*2+1])+1./(geo.w[(i+1)*2+1]*geo.w[(i+1)*2+1]));
        
   geo.w[2*configs.freq_num_fd-1] = 0.001;
    
//Read td noise
    
   fgets(buf, 2000, conf);
   fgets(buf, 2000, conf);
   parse(buf, values);
    
   for(i = 0; i < configs.num_channels; i++)
   {
      geo.w[2*configs.freq_num_fd+i] = 1/values[i];
      geo1.w[i] = 1/values[i];
   }
   
}
