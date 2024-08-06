#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <stdlib.h>
#include <string.h>

#include "fft23.h"

#define SQRT3D2 .8660254037844386

void inv_ind(FFT *fft) 
{
    int i,j,k,l,m,num = 1,offs,order = 0;
    int deo = !(fft->order2);
    for(i=0;i<fft->n;i++)
	fft->inv_index[i] = i;

    for(;order<fft->order3-deo;order++) {
	l = (order&1);
	k = !l;
	if(k) k*=fft->n;
	else  l*=fft->n;
	m = fft->n/num - 1;
	for(i=0;i<num;i++) {
	    offs = i*(m+1);
	    for(j=0;j<m;j++) {
		fft->inv_index[j+offs+k] = fft->inv_index[(j*3)%m+offs+l];
	    }
	    fft->inv_index[m+offs+k] = fft->inv_index[m+offs+l];
	}
	num *= 3;
    }

 
    for(;order<fft->order2+fft->order3-1;order++) {
	l = (order&1);
	k = !l;
	if(k) k*=fft->n;
	else  l*=fft->n;
	m = fft->n/num - 1;
	for(i=0;i<num;i++) {
	    offs = i*(m+1);
	    for(j=0;j<m;j++) {
		fft->inv_index[j+offs+k] = fft->inv_index[(j*2)%m+offs+l];
	    }
	    fft->inv_index[m+offs+k] = fft->inv_index[m+offs+l];
	}
	num *= 2;
    }
    if(k) 
	memcpy(fft->inv_index,&fft->inv_index[k],k*sizeof(int));
}

void init_fft(FFT *fft, int n)
{
    int nn = n,order = 1;
    fft->order2 = fft->order3 = 0;
    fft->n = 1;

    while (nn&&(!(nn&1))) {
	nn >>= 1;
	fft->order2++;
	fft->n <<= 1;
    }

    while (nn&&(!(nn%3))) {
	nn/=3;
	fft->order3++;
	fft->n *= 3;
    }

    fft->inv_index = (int *)realloc(fft->inv_index,2*fft->n*sizeof(int));
    fft->wkn= (double complex *)realloc(fft->wkn,fft->n*sizeof(double complex));
    fft->xn = (double complex *)realloc(fft->xn,fft->n*sizeof(double complex));
    fft->fn = (double complex *)realloc(fft->fn,fft->n*sizeof(double complex));
    
    memset(fft->inv_index,0,2*fft->n*sizeof(int));
    memset(fft->wkn,0,fft->n*sizeof(double complex));
    memset(fft->xn,0,fft->n*sizeof(double complex));
    memset(fft->fn,0,fft->n*sizeof(double complex));

    int i,j;

    inv_ind(fft);

    for(i=0;i<fft->n;i++) {
	fft->wkn[i] = cexp(-2.*i*I*M_PI/fft->n);
    }
}

int wkn(int k, int n_c, FFT *fft)
{
return (k*(fft->n/n_c));
}

FFT* FFT_new(int w)
{
    FFT *f;
    f=(FFT*)malloc(sizeof(FFT));
    memset(f,0,sizeof(FFT));
    init_fft(f,w);
    return f;
}

void FFT_free(FFT *f)
{
    free(f);
}

void FFT_free_full(FFT *f)
{
    free(f->inv_index);
    free(f->wkn);
    free(f->xn);
    free(f->fn);
    //free(f);
}

void fft_pro(FFT *fft,int inv)
{
    int n, nd2, k, mpnd2, mpnd3, m, order;
    double complex temp,temp1,temp2;

for(m=0;m<fft->n;m++)
	fft->fn[m] = fft->xn[fft->inv_index[m]];

for(n=1,order=0;order<fft->order2;order++) {
    nd2 = n, n+= n;
    for(k=0;k<nd2;k++) {
	int ind = wkn(k,n,fft);
	for (m=k; m<fft->n; m+=n) {
	    mpnd2 = m + nd2;
        temp = fft->wkn[ind] * fft->fn[mpnd2];
        fft->fn[mpnd2] = fft->fn[m] - temp;
	    fft->fn[m]     = fft->fn[m] + temp;
	}
    }
    
}

for(;order<fft->order2+fft->order3;order++) {
    nd2 = n, n+= 2*n;
    for(k=0;k<nd2;k++) {
	int ind  = wkn(k,n,fft);
	int ind1 = wkn(2*k,n,fft);
	for (m=k; m<fft->n; m+=n) {
	    mpnd2 = m + 2*nd2;
	    mpnd3 = m + nd2;
	    temp  = fft->wkn[ind] * fft->fn[mpnd3];
	    temp1 = fft->wkn[ind1]* fft->fn[mpnd2];
	    temp2 = temp + temp1;   // sum
	    temp  = temp - temp1;  
	    temp1 = I*(creal(temp)*SQRT3D2) -   (cimag(temp)*SQRT3D2);  // i*diff*sqrt3/2
	    temp  =   (.5*creal(temp2))     + I*(.5*cimag(temp2))    ;  // half sum
	    fft->fn[mpnd3] = fft->fn[m] - temp - temp1;
	    fft->fn[mpnd2] = fft->fn[m] - temp + temp1;
	    fft->fn[m]     = fft->fn[m] + temp2;
	}
    }
}

if(!inv) {
    for(m=1;m<fft->n/2;m++)
        {
        fft->fn[m] += conj(fft->fn[fft->n - m]);
        fft->fn[fft->n - m] = 0;
        fft->fn[m] = conj(fft->fn[m])/fft->n;
        }
    fft->fn[0] = 0;
    fft->fn[fft->n/2] = 0;
} else {
    ;//for(m=0;m<fft->n;m++)
     //   fft->fn[m] /= fft->n;
}

}

/*

FFT fft,ifft;

int main()
{
#define fre (200000./2592.*2048.)
#define dim 2592
//2592    
// 81 * 16 !!! 81 * 32 &&&???
#define len 64
//81

    int i,ln;
    double j,dj;

    double freq;
    double waveform[3000];
    char buf[100];
    
    FILE *fout, *fin;
    
    fin = fopen("waveform.XYZ", "rt");
    
    for(i=0;fgets(buf,100,fin);i++){
        waveform[i] = atof(buf);
    }
    
    fout = fopen("res.XYZ", "wt");

    printf("/ ind  sign inv_sign  freq spec re im \n");
    fprintf(fout, "ind  sign inv_sign  freq spec re im \n");

    memset(&fft,0,sizeof(fft));
    memset(&ifft,0,sizeof(ifft));
    init_fft(&fft,dim);
    init_fft(&ifft,dim);
    for(i=0,j=0;i<dim;i++) {
        fft.xn[i] = waveform[i];
    }

    fft_pro(&fft,0);
    for(i=0;i<ifft.n;i++) {
	ifft.xn[i] = fft.fn[i];
    }
    fft_pro(&ifft,1);
    double arg, argp = 0, sign = 1.;
    for(i=0;i<dim;i+=1) {
	freq = fre/dim*i;
	arg = carg(ifft.fn[i])/M_PI*180.;
	argp = arg;
	//printf("%03d  %lf %lf  %lf  %lf %lf %lf  \n ",i,creal(fft.xn[i]),creal(ifft.fn[i]), 
	//       freq,cabs(fft.fn[i])/fft.n,fft.fn[i]/fft.n);
	       
	//fprintf(fout, "%03d  %lf %lf  %lf  %lf %lf %lf  \n ", i, creal(fft.xn[i]), creal(ifft.fn[i]), 
	//       freq, cabs(fft.fn[i])/fft.n, fft.fn[i]/fft.n);
	fprintf(fout, "%d %lf %lf %lf %lf %lf %lf \n ", ifft.inv_index[i], ifft.wkn[i], ifft.xn[i], ifft.fn[i]);
    }
    fclose(fout);
return 0;
}
*/

