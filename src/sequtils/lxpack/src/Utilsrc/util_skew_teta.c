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


#define WINDOW_DFT  1000
#define STEP_DFT    300

/* ----------------------------------------------- */
static int sCount(char *seq, int from, int to, int len, int c) {

  int i, count = 0;
  for (i = from ; i < to ; i++) {
    if (seq[i%len] == c) 
      count++;
  } 
  return count;
}

/* ----------------------------------------------- */
/*ARGSUSED*/

main(argn, argv)
    int  argn;
    char *argv[];
{
    int nbseq, opt, len, pos;
    int ng[2], nc[2], dg[2], dc[2];
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
        
        /* rotate origin */

        for (pos = 0 ; pos < len; pos += step) {

          if ((pos == 0) || (width < step)) {
        
            ng[0] = sCount(seq->seq, pos, pos + width, len, 'G');
            nc[0] = sCount(seq->seq, pos, pos + width, len, 'C');
        
            ng[1] = sCount(seq->seq, pos + len - width, pos + len, len, 'G');
            nc[1] = sCount(seq->seq, pos + len - width, pos + len, len, 'C');

          }
          else {
          
            dg[0] = sCount(seq->seq, pos - step, pos, len, 'G');
            dc[0] = sCount(seq->seq, pos - step, pos, len, 'C');
            
            dg[1] = sCount(seq->seq, pos - step + width, pos + width, len, 'G');
            dc[1] = sCount(seq->seq, pos - step + width, pos + width, len, 'C');
            
            ng[0] += dg[1] - dg[0];
            nc[0] += dc[1] - dc[0];
            
            dg[1] = sCount(seq->seq, pos - step + len - width, pos + len - width, len, 'G');
            dc[1] = sCount(seq->seq, pos - step + len - width, pos + len - width, len, 'C');
          
            ng[1] += dg[0] - dg[1];
            nc[1] += dc[0] - dc[1];
          }
        
          skew =   ((float) (ng[0] - nc[0]) / (float) (ng[0] + nc[0]))
                 - ((float) (ng[1] - nc[1]) / (float) (ng[1] + nc[1]));
            
          printf("%d %f\n", pos, skew);
       }
    }

    FreeFastaSequence(seq);

    exit(0);
    
    /*NOTREACHED*/
}
