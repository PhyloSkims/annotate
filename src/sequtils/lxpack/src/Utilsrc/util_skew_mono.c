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
#define STEP_DFT    10000

#define NB_NUC 4

static char sDna[] = "ACGT";

/* ----------------------------------------------- */
static char *sIndexSeq(char *seq, char *alpha) {

  int len;
  char *c, *pos;
  
  len = strlen(alpha);
  
  for (c = seq; *c ; c++) {
    if (pos = strchr(alpha, *c))
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
/*ARGSUSED*/

main(argn, argv)
    int  argn;
    char *argv[];
{
    int opt, len, pos, nuc;
    int countLeft[NB_NUC + 1], countRight[NB_NUC + 1];
    int window, width, step;
    float skew;
    FastaSequence *seq;

    extern char *optarg;    /* externs for getopts (3C) */

    step   = STEP_DFT;      /* default step size   */
    window = WINDOW_DFT;    /* defaut window = 50% */
    
    while ((opt = getopt(argn, argv, "hs:w:W:")) != -1) {

        switch (opt) {

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
    
    seq = NewFastaSequence();

    ReadFastaSequence(stdin, seq);

		if (! seq->ok) {
		  (void) printf("error while reading sequence\n");
		  exit(3);
		}

		len = seq->length;

		if (window > 0)
			width = window;
		else
			width = len * (-window) / 100;
		
		/* index sequence */
		
		(void) sIndexSeq(seq->seq, sDna);
		
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
		
      printf("%d", pos);
      for (nuc = 0 ; nuc < NB_NUC ; nuc++) {
        skew = (float) (countRight[nuc] + countLeft[nuc]);
        if (skew == 0.) 
          skew = 1.;
        else
          skew =  (float) (countRight[nuc] - countLeft[nuc]) / skew;
        printf(" %f6.2", skew * 100.);
      }
			printf("\n");
			
	 }

   FreeFastaSequence(seq);

   exit(0);
}
