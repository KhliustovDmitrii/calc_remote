/****************************\
 * Uncertainty measures     *
 * Currently implemented:   *
 * variance based measure   *
 * resistivity based measure*
\****************************/

//Variance based measure
void var_mes(double *S0, double *S1, double *mes_buf, int sdim)
{
   int i, j;
   double v, v1;
   for(i=0;i<sdim;i++) 
   {
      v = 0;
      v1 = 0;
      for (j=i;j<sdim;j++) 
      {
         v+= S0[i*sdim+j]*S0[i*sdim+j];
         v1+= S1[i*sdim+j]*S1[i*sdim+j];
      }
      mes_buf[i] = 1 - sqrt(v/v1);
   }
}

//Resistivity based measure
void rho_mes(double rho_ini, double *rho, double *mes_buf, int lay_num)
{
   int i;
   for(i=0; i<lay_num; i++)
      mes_buf[i] = (2/M_PI)*atan(fabs(log(rho[i]) - log(rho_ini))/fabs(log(rho[i]) + log(rho_ini)));
}
