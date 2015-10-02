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

#ifdef SGI
#include <getopt.h>
#endif

#include "libbio.h"
#include "libfasta.h"


#define WINDOW_DFT  -25 /* aka 25% */
#define STEP_DFT    10000

#define NB_NUC 4

#define A 0
#define C 1
#define G 2
#define T 3

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
    int nbseq, opt, len, pos;
    int window, width, step;
    int count[NB_NUC];
    float skew;
    FastaSequence *seq;

    extern char *optarg;    /* externs for getopts (3C) */

    window = WINDOW_DFT;    /* default window size */
    step   = STEP_DFT;      /* default step size   */
    
    while ((opt = getopt(argn, argv, "hs:w:W:")) != -1) {

        switch (opt) {

        case 'h':
            (void) printf("compute GC skew\n");
            (void) printf("usage: util_skew [-w width] [-s step]\n");
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
            (void) printf("usage: util_skew [-w width] [-s step]\n");
            exit(6);
            break;
        }
    }
    
    seq = NewFastaSequence();

    nbseq = 0;

    while (ReadFastaSequence(stdin, seq)) {

        nbseq++;

        if (! seq->ok) {
        (void) printf("error at seq # %d\n", nbseq);
        continue;
        }

        len = seq->length - width;

				if (window > 0)
					width = window;
				else
					width = len * (-window) / 100;

		    (void) sIndexSeq(seq->seq, sDna);
        
        skew = 0.;
        
        for (pos = 0 ; pos < len ; pos += step) {

					if ((pos == 0) || (width < step)) {
						sCount(seq->seq, pos + len - width, pos + len + width, len, count, 0);
					}
					else {
						sCount(seq->seq, pos - step + len - width, pos + len - width, len, count, -1);
						sCount(seq->seq, pos - step + width, pos + width, len, count, 1);
					}

          skew += (float) (count[G] - count[C]) / (float) (count[G] + count[C]);

          printf("%d %f\n", pos, skew);
          
        }

    }

    FreeFastaSequence(seq);

    exit(0);
    
    /*NOTREACHED*/
}
