/* ---------------------------------------------------------------- */
/* Copyright (c) Atelier de BioInformatique                         */
/* @file: util_skew.c                                               */
/* @desc: compute GC skew                                           */
/*                                                                  */
/* @history:                                                        */
/* @+       <Gloup> : Jan 96 : PWG version                          */
/* ---------------------------------------------------------------- */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifdef SGI
#include <getopt.h>
#endif

#include "libbio.h"
#include "libfasta.h"

#define WINDOW_DFT  -25  /* aka 25% */
#define STEP_DFT    300

#define NB_NUC 4

static char sDna[] = "ACGT";

/* ----------------------------------------------- */
/* remove invalid and duplicate symbols in         */
/* alphabet                                        */
static int sCleanAlpha(char *alpha, char *clean) {
  char *c, *cc;
  
  for (c = alpha, cc = clean; *c ; c++) {
    if (    (strchr(sDna, *c) != NULL)
         && (strchr(alpha, *c) == c))
      *cc++ = *c;
  }
  *cc = '\000';
  return strlen(clean);
}


/* ----------------------------------------------- */
static char *sIndexSeq(char *seq, char *alpha) {

  int len;
  char *c, *pos;
  
  len = strlen(alpha);
  
  for (c = seq; *c ; c++) {
    if (strchr(sDna, *c) == NULL)
      *c = -1;
    else if ((pos = strchr(alpha, *c)) != NULL)
      *c = (int) (pos - alpha);
    else
      *c = len;
  }
  return seq;
}

/* ----------------------------------------------- */
static void sCount(char *seq, int from, int to, int len, int *count, int factor) {

  int i, j;
  
  if (factor == 0) {
    for (i = 0 ; i < NB_NUC ; i++)
      count[i] = 0;
    factor = 1;
  }
  
  for (i = from ; i < to ; i++) {
    j = seq[i%len]; 
    if (j >= 0) 
      count[j] += factor;
  } 
}

/* ----------------------------------------------- */
/* ----------------------------------------------- */

#define SQRTF(x) (float) sqrt((float) x)
#define FABSF(x) (float) fabs((float) x)


/* ---------------------- */

/* Numerical Recipes standard error handler */

static void nrerror(char *error_text)
{
  fprintf(stderr, "Numerical Recipes run-time error...\n");
  fprintf(stderr, "%s\n", error_text);
  exit(1);
}

/* ---------------------- */

static float gammln(float xx)
{
  double x, y, tmp, ser;
  static double cof[6] = {
    76.18009172947146, -86.50532032941677,
    24.01409824083091, -1.231739572450155,
    0.1208650973866179e-2, -0.5395239384953e-5
  };
  int j;

  y = x = xx;
  tmp = x + 5.5;
  tmp -= (x+0.5) * log(tmp);
  ser = 1.000000000190015;
  for (j = 0 ; j <= 5 ; j++)
    ser += cof[j] / ++y;
    
  return -tmp + log(2.5066282746310005 * ser / x);
}

/* ---------------------- */

#define ITMAX 100
#define EPS 3.0e-7
#define FPMIN 1.0e-30

static void gcf(float *gammcf, float a, float x, float *gln)
{
  int i;
  float an, b, c, d, del, h;

  *gln = gammln(a);
  b = x + 1.0 - a;
  c = 1.0/FPMIN;
  d = 1.0/b;
  h = d;
  for (i = 1 ; i <= ITMAX ; i++) {
    an = -i*(i-a);
    b += 2.0;
    d = an*d + b;
    if (FABSF(d) < FPMIN)
      d = FPMIN;
    c = b + an/c;
    if (FABSF(c) < FPMIN)
      c=FPMIN;
    d = 1.0/d;
    del = d*c;
    h *= del;
    if (FABSF(del-1.0) < EPS)
      break;
  }
  if (i > ITMAX)
    nrerror("a too large, ITMAX too small in gcf");
    
  *gammcf = exp (-x+a*log(x)-(*gln)) * h;
}
#undef ITMAX
#undef EPS
#undef FPMIN

/* ---------------------- */

#define ITMAX 100
#define EPS 3.0e-7

static void gser(float *gamser, float a, float x, float *gln)
{
  int n;
  float sum, del, ap;

  *gln = gammln(a);
  if (x <= 0.0) {
    if (x < 0.0)
      nrerror("x less than 0 in routine gser");
    *gamser=0.0;
    return;
  } else {
    ap = a;
    del = sum = 1.0/a;
    for (n = 1 ; n <= ITMAX ; n++) {
      ++ap;
      del *= x/ap;
      sum += del;
      if (FABSF(del) < FABSF(sum)*EPS) {
        *gamser = sum * exp(-x+a*log(x)-(*gln));
        return;
      }
    }
    nrerror("a too large, ITMAX too small in routine gser");
    return;
  }
}
#undef ITMAX
#undef EPS

/* ---------------------- */

static float gammq(float a, float x)
{
  float gamser,gammcf,gln;

  if (x < 0.0 || a <= 0.0)
    nrerror("Invalid arguments in routine gammq");
  if (x < (a+1.0)) {
    gser(&gamser, a, x, &gln);
    return 1.0-gamser;
  } else {
    gcf(&gammcf, a, x, &gln);
    return gammcf;
  }
}

/* ---------------------- */

void chstwo(int bins1[], int bins2[], int nbins, 
            float *df, float *chsq, float *prob)
{
  int j;
  float temp, s1, s2, c1, c2;

  s1 = s2 = 0.;
  for (j=0 ; j<nbins ; j++) {
    s1 += bins1[j];
    s2 += bins2[j];
  }
    
  c1 = SQRTF(s2/s1);
  c2 = SQRTF(s1/s2);

  *df = nbins-1;
  *chsq = 0.0;
  for (j = 0 ; j < nbins ; j++) {
    if ((bins1[j] == 0) && (bins2[j] == 0))
      --(*df);
    else {
      temp = bins1[j]*c1 - bins2[j]*c2;
      *chsq += temp*temp/(bins1[j]+bins2[j]);
    }
  }
  *prob = gammq(0.5*(*df),0.5*(*chsq));
}

/* ----------------------------------------------- */
/*ARGSUSED*/

main(argn, argv)
    int  argn;
    char *argv[];
{
    int nbseq, opt, len, nbins, pos, sign;
    int countLeft[NB_NUC], countRight[NB_NUC];
    int window, width, step;
    char alpha[NB_NUC + 1];
    float chi, proba, df;
    FastaSequence *seq;

    extern char *optarg;    /* externs for getopts (3C) */

    step   = STEP_DFT;      /* default step size   */
    window = WINDOW_DFT;    /* defaut window = 50% */
    
    (void) strcpy(alpha, sDna);
    
    while ((opt = getopt(argn, argv, "a:hs:w:W:")) != -1) {

        switch (opt) {

		    case 'a':
		        if (   (strlen(optarg) > NB_NUC)
		            || (sCleanAlpha(optarg, alpha) == 0)) {
		          (void) printf("bad alpha value\n");
               exit(5);
            }
		        break;

        case 'h':
            (void) printf("compute GC skew\n");
            (void) printf("usage: util_skew [-w|W width] [-s step]\n");
            exit(0);
            break;

        case 's':
            if (   (sscanf(optarg, "%d", &step) != 1)
            || (step <= 0) ) {
               (void) printf("bad step value\n");
               exit(5);
            }
            break;

        case 'w':
            if (   (sscanf(optarg, "%d", &window) != 1)
            || (window <= 0) ) {
               (void) printf("bad window value\n");
               exit(5);
            }
            break;

        case 'W':
            if ( (sscanf(optarg, "%d", &window) != 1)
            || (window <= 0) 
            || (window > 100)) {
               (void) printf("bad window percent value\n");
               exit(5);
            }
            window = -window;
            break;

        case '?':
            (void) printf("usage: util_skew [-w|W width] [-s step]\n");
            exit(6);
            break;
        }
    }
    
    nbins = strlen(alpha);
    if (nbins < NB_NUC)
      nbins++;
    
    seq = NewFastaSequence();

    nbseq = 0;

    while (ReadFastaSequence(stdin, seq)) {

        nbseq++;

        if (! seq->ok) {
        (void) printf("error at seq # %d\n", nbseq);
        continue;
        }

        len = seq->length;

        if (window > 0)
          width = window;
        else
          width = len * (-window) / 100;
        
        /* index sequence */
        
        (void) sIndexSeq(seq->seq, alpha);
        
        /* rotate origin */

        for (pos = 0 ; pos < len; pos += step) {

          if ((pos == 0) || (width < step)) {
        
            sCount(seq->seq, pos, pos + width, len, countRight, 0);
            sCount(seq->seq, pos + len - width, pos + len, len, countLeft, 0);
          }
          else {
          
            sCount(seq->seq, pos - step, pos, len, countRight, -1);
            sCount(seq->seq, pos - step + width, pos + width, len, countRight, 1);
            
            sCount(seq->seq, pos - step + len - width, pos + len - width, len, countLeft, -1);
            sCount(seq->seq, pos - step, pos, len, countLeft, 1);
          }
        
          chstwo(countRight, countLeft, nbins, &df, &chi, &proba);

          /* add a sign */
          
          sign = ((nbins == 2) ? ((countRight[0] >= countLeft[0]) ? 1 : -1)
                              : 1);
          printf("%d %f %g\n", pos, chi * sign, proba);
       }
    }

    FreeFastaSequence(seq);

    exit(0);
    
    /*NOTREACHED*/
}
