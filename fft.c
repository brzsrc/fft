/* Factored discrete Fourier transform, or FFT, and its inverse iFFT */
/* from https://www.math.wustl.edu/~victor/mfmm/fourier/fft.c */

/* define TEST to print test results  */
/* #define TEST */

/* Compile with gcc -DSCALE=16 for 16-point transform. */

/* define OURSINCOS to use our own implementation of sin and cos (slightly approximate) */

/* eg gcc -DOURSINCOS -DSCALE=16 fft.c sincos.c */
/* or gcc -DSCALE=16 fft.c -lm to use the library sin and cos */

/* To compile for simplescalar you need to add -DNOFPCONSTANTS because simplescalar's gcc can't handle 
 * floating point constants (!).  For simplescalar: 
 * gcc -DOURSINCOS -DNOFPCONSTANTS -DSCALE=16 fft.c sincos.c
 */

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#ifndef SCALE
#define SCALE	3		/* for 2^3 points */
#endif

#ifdef NOFPCONSTANTS
int N;
float fN;
#else
#define N	(1<<SCALE)		/* N-point FFT, iFFT */
#define fn (float)(N)
#endif

typedef float real;
typedef struct{real Re; real Im;} complex;

#ifdef NOFPCONSTANTS
extern float PI; /* actually dclared in sincos.c */
float FLT_MIN;
float f2;
#else
/* #define PI   3.14159265358979323846264338327950288 */
#define PI      3.14159265358
#ifndef FLT_MIN
#define FLT_MIN -3.40282346638528859811704183484516925440e+38
#define f2 2.0
#endif
#endif




/* for simplescalar we need our own implementation of sin and cos */
#ifdef OURSINCOS
float thecos(float);
float thesin(float);
#else
#define thecos cos
#define thesin sin
#endif

/* Print a vector of complexes as ordered pairs. */
static void
print_vector(
	     const char *title,
	     complex *x,
	     int n)
{
  int i;
  printf("%s (dim=%d):", title, n);
  for(i=0; i<n; i++ ) printf(" %5.2f,%5.2f ", x[i].Re,x[i].Im);
  putchar('\n');
  return;
}

/* Find index of maximum real value in vector. */
static int
Re_imax_vector(
	     complex *x,
	     int n)
{
  int i;
  int imax = -1;
  float max = FLT_MIN;
  for(i=0; i<n; i++ ) {
    if (x[i].Re > max) {
      max = x[i].Re;
      imax = i;
    }
  }
  return imax;
}

/* Find maximum real value in vector. */
static float
Re_max_vector(
	     complex *x,
	     int n)
{
  int i;
  float max = FLT_MIN;
  for(i=0; i<n; i++ ) {
    if (x[i].Re > max) {
      max = x[i].Re;
    }
  }
  return max;
}


/* 
   fft(v,N):
   [0] If N==1 then return.
   [1] For k = 0 to N/2-1, let ve[k] = v[2*k]
   [2] Compute fft(ve, N/2);
   [3] For k = 0 to N/2-1, let vo[k] = v[2*k+1]
   [4] Compute fft(vo, N/2);
   [5] For m = 0 to N/2-1, do [6] through [9]
   [6]   Let w.re = cos(2*PI*m/N)
   [7]   Let w.im = -sin(2*PI*m/N)
   [8]   Let v[m] = ve[m] + w*vo[m]
   [9]   Let v[m+N/2] = ve[m] - w*vo[m]
 */
void
fft( complex *v, int n, complex *tmp )
{
  if(n>1) {			/* otherwise, do nothing and return */
    int k,m;    complex z, w, *vo, *ve;
    ve = tmp; vo = tmp+n/2;
    for(k=0; k<n/2; k++) {
      ve[k] = v[2*k];
      vo[k] = v[2*k+1];
    }
    fft( ve, n/2, v );		/* FFT on even-indexed elements of v[] */
    fft( vo, n/2, v );		/* FFT on odd-indexed elements of v[] */
    for(m=0; m<n/2; m++) {
      w.Re = thecos(f2*PI*m/(double)n);
      w.Im = -thesin(f2*PI*m/(double)n);
      z.Re = w.Re*vo[m].Re - w.Im*vo[m].Im;	/* Re(w*vo[m]) */
      z.Im = w.Re*vo[m].Im + w.Im*vo[m].Re;	/* Im(w*vo[m]) */
      v[  m  ].Re = ve[m].Re + z.Re;
      v[  m  ].Im = ve[m].Im + z.Im;
      v[m+n/2].Re = ve[m].Re - z.Re;
      v[m+n/2].Im = ve[m].Im - z.Im;
    }
  }
  return;
}

/* 
   ifft(v,N):
   [0] If N==1 then return.
   [1] For k = 0 to N/2-1, let ve[k] = v[2*k]
   [2] Compute ifft(ve, N/2);
   [3] For k = 0 to N/2-1, let vo[k] = v[2*k+1]
   [4] Compute ifft(vo, N/2);
   [5] For m = 0 to N/2-1, do [6] through [9]
   [6]   Let w.re = cos(2*PI*m/N)
   [7]   Let w.im = sin(2*PI*m/N)
   [8]   Let v[m] = ve[m] + w*vo[m]
   [9]   Let v[m+N/2] = ve[m] - w*vo[m]
 */
void
ifft( complex *v, int n, complex *tmp )
{
  if(n>1) {			/* otherwise, do nothing and return */
    int k,m;    complex z, w, *vo, *ve;
    ve = tmp; vo = tmp+n/2;
    for(k=0; k<n/2; k++) {
      ve[k] = v[2*k];
      vo[k] = v[2*k+1];
    }
    ifft( ve, n/2, v );		/* FFT on even-indexed elements of v[] */
    ifft( vo, n/2, v );		/* FFT on odd-indexed elements of v[] */
    for(m=0; m<n/2; m++) {
      w.Re = thecos(f2*PI*m/(double)n);
      w.Im = thesin(f2*PI*m/(double)n);
      z.Re = w.Re*vo[m].Re - w.Im*vo[m].Im;	/* Re(w*vo[m]) */
      z.Im = w.Re*vo[m].Im + w.Im*vo[m].Re;	/* Im(w*vo[m]) */
      v[  m  ].Re = ve[m].Re + z.Re;
      v[  m  ].Im = ve[m].Im + z.Im;
      v[m+n/2].Re = ve[m].Re - z.Re;
      v[m+n/2].Im = ve[m].Im - z.Im;
    }
  }
  return;
}


int
main(void)
{
  /*   complex v[N], v1[N], v2[N], scratch[N]; */
  complex *v = malloc(N*sizeof(complex));
  complex *v1 = malloc(N*sizeof(complex));
  complex *v2 = malloc(N*sizeof(complex));
  complex *scratch = malloc(N*sizeof(complex));
  int k;
#ifdef NOFPCONSTANTS
  int n2=2, n3=3, n4=4, n8=8, n29=29, n44=44, n60=60;
  float f4=(float)(n4);
  /* base 60 turns out to be good for pi, see wikipedia */
  PI=(float)(n3)+(float)(n8)/(float)(n60)+(float)(n29)/(float)(n60*n60)+(float)(n44)/(float)(n60*n60*n60);
  f2 = (float)(n2);
  N = (1<<SCALE); fN = (float)(N);
  FLT_MIN = f2/(f2*f2*f2*f2*f2*f2*f2*f2*f2*f2); FLT_MIN *= FLT_MIN; /* smallish */
  /* printf("pi=%f f2=%f flt_min=%f\n", PI, f2, FLT_MIN); */
#endif
  printf("Size = %d\n", N);

  /*
  printf("sin(PI/2=%f) = %f\n", PI/f2, thesin(PI/f2));
  printf("cos(PI/2=%f) = %f\n", PI/f2, thecos(PI/f2));
  printf("sin(PI/4=%f) = %f\n", PI/f4, thesin(PI/f4));
  printf("cos(PI/4=%f) = %f\n", PI/f4, thecos(PI/f4));
  printf("cos(0) = %f\n", thecos(PI-PI));
  exit(0);
  */
  
  /* Fill v[] with a function of known FFT: */
  for(k=0; k<N; k++) {
#ifdef NOFPCONSTANTS
    int n1=1, n8=8, n3=3, n10=10;
    float f1=(float)(n1), f8=(float)(n8), f3=(float)(n3), f10=(float)(n10), fpoint3=f3/f10, feighth=f1/f8;
#else
    #define feighth 0.125
    #define fpoint3 0.3
    #define f2 2
#endif
    v[k].Re = feighth*thecos(f2*PI*k/fN);
    /*    printf("feighth=%f fpoint3=%f f2=%f f2*PI*k/fN = %f, cos()=%f\n", feighth, fpoint3, f2, f2*PI*k/fN, thecos(f2*PI*k/fN));*/
    v[k].Im = feighth*thesin(f2*PI*k/fN);
#ifdef FALSE
    v1[k].Re =  fpoint3*thecos(f2*PI*k/fN);
    v1[k].Im = -fpoint3*thesin(f2*PI*k/fN);

    v2[k].Re =  fpoint3*(thecos(f2*PI*k/fN) + fpoint3*thecos(f3*PI*k/fN));
    v2[k].Im = -fpoint3*(thesin(f2*PI*k/fN) - fpoint3*thesin(f3*PI*k/fN));
#endif
  }

#ifdef TEST
  /* FFT, iFFT of v[]: */
  print_vector("Orig", v, N);
  fft( v, N, scratch );
  print_vector(" FFT", v, N);
  ifft( v, N, scratch );
  print_vector("iFFT", v, N);

  /* FFT, iFFT of v1[]: */
  /*  print_vector("Orig", v1, N);
  fft( v1, N, scratch );
  print_vector(" FFT", v1, N);
  ifft( v1, N, scratch );
  print_vector("iFFT", v1, N); */

  /* FFT, iFFT of v1[]: */
  /*  print_vector("Orig", v2, N);
  fft( v2, N, scratch );
  print_vector(" FFT", v2, N);
  ifft( v2, N, scratch );
  print_vector("iFFT", v2, N); */
#else

  printf("Max            v[i].Re is at %d (%f)\n", Re_imax_vector( v, N ), Re_max_vector( v, N ));
  fft( v, N, scratch );
  printf("Max       fft(v)[i].Re is at %d (%f)\n", Re_imax_vector( v, N ), Re_max_vector( v, N ));
  ifft( v, N, scratch );
  printf("Max ifft(fft(v))[i].Re is at %d (%f)\n", Re_imax_vector( v, N ), Re_max_vector( v, N ));

#endif

  
  /*  exit(EXIT_SUCCESS); */
}
