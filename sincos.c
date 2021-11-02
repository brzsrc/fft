/* from https://web.eecs.utk.edu/~azh/blog/cosine.html (cos_taylor_running_yterms) */

/* This is a butchered a bit to eliminate all floating-point constants, because the simplescalar compiler can't handle them.  
 * This is activated by the macro NOFPCONSTANTS.
 */

#ifdef SINCOSTEST
#include <math.h>
#include <stdio.h>
#define thecos cos
#define thesin sin
#else
#define approxcos thecos
#define approxsin thesin
#endif

#define Y 8
#define TYPE float

#ifdef NOFPCONSTANTS
TYPE PI;
#else
/* #define PI	3.14159265358979323846264338327950288 */
#define PI	3.14159265358 
#endif

TYPE approxcos(TYPE x)
{
    int div = (int)(x / PI);
    char sign = 1;
#ifdef NOFPCONSTANTS
    int n1 = 1;
    TYPE f1 = (TYPE)(n1);
    TYPE result = f1;
    TYPE inter = f1;
#else
    TYPE result = 1.0;
    TYPE inter = 1.0;
#endif
    TYPE num = x * x;
    int i;

    x = x - (div * PI);
    if (div % 2 != 0)
        sign = -1;

    for (i = 1; i <= Y; i++)
    {
#ifdef NOFPCONSTANTS
      TYPE comp = (float)(i+i);
      TYPE den = comp * (comp - f1);
#else
      TYPE comp = 2.0 * i;
      TYPE den = comp * (comp - 1.0);
#endif
      /*      printf("i=%d comp = %f  den = %f\n", i, comp, den);*/
        inter *= num / den;
        if (i % 2 == 0)
            result += inter;
        else
            result -= inter;
    }
    /*    printf("approxcos(%f)=%f\n", x, sign*result); */
    
    return sign * result;
}

TYPE approxsin(TYPE x) { int n2=2; TYPE f2=(float)(n2); return approxcos(PI/f2 - x); }


#ifdef SINCOSTEST

/* #define DELTA (PI/20.0f) */
#ifdef NOFPCONSTANTS
TYPE DELTA;
#else
#define DELTA (TYPE)0.314159
#endif

int main(){
  TYPE x;
  int i;

#ifdef NOFPCONSTANTS
  int n0 = 0, n10=10;
  int n2 = 2; TYPE f2 = (float)(n2);
  int n4 = 2; TYPE f4 = (float)(n4);
  int n22 = 22;
  int n7 = 7;
  /*   PI=(float)(n22)/(float)(n7); */
  int n3=3, n8=8, n29=29, n44=44, n60=60;
  PI=(float)(n3)+(float)(n8)/(float)(n60)+(float)(n29)/(float)(n60*n60)+(float)(n44)/(float)(n60*n60*n60); /* base 60 turns out to be good for pi, see wikipedia */
  DELTA = PI/(float)(n10);
  printf("PI= %f DELTA = %f\n", PI, DELTA);
#endif

  printf("sin(PI/2) = %f\n", approxsin(PI/f2));
  printf("cos(PI/2) = %f\n", approxcos(PI/f2));
  printf("sin(PI/4) = %f\n", approxsin(PI/f4));
  printf("cos(PI/4) = %f\n", approxcos(PI/f4));
  printf("cos(0) = %f\n", approxcos(PI-PI));

#ifdef NOFPCONSTANTS
  x = (float)n0;
#else
  x = 0.0;
#endif
  
#ifdef FALSE
  for (i=0; i<20; ++i) {
    x = x+DELTA;
    printf("%d %f\n", i, x);
  }
#endif

  for (x = -PI; x < PI; x+=DELTA) {
    printf("cos(%f) = %f, sin(%f) = %f\n", x, approxcos(x), x, approxsin(x));
    printf("%f\n", x);
    /*    printf("cos(%f) = %f (vs %f), sin(%f) = %f (vs %f)\n", x, approxcos(x), cos(x), x, approxsin(x), sin(x));*/
  }
}
#endif
