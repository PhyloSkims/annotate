/* ==================================================== */
/* @file: kimfit_main.c                                 */
/* triangle fit program                                 */
/* @history:                                            */
/* @+   Apr. 97 <AV> first draft                        */
/* ==================================================== */

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "kimono.h"

#ifndef BUFSIZ
#define BUFSIZ          1024
#endif

#define NBEST 2


/* ----------------------------------------------- */
/* look for best window on signal                  */
/* ----------------------------------------------- */

static float sLookBest(float *x, int n, int iwidth) {
  int i;
  float sum, best;
  
  if (iwidth > n) iwidth = n;
  
  sum = 0.;
  
  for (i = 0 ; i < iwidth ; i++) {
    sum += x[i];
  }

  best = sum;

  for (i = 1 ; i < n ; i++) {

    sum -= x[i-1];
    sum += x[(i+iwidth-1) % n];

    if (sum > best)
      best = sum;

  }
  
  return best;
}

/* ----------------------------------------------- */
/* fit triangle on signal                          */
/* ----------------------------------------------- */

static int sFitTriangle(float *x, int n, int iwidth, float height) {
  int i, j, fitpos;

  if (iwidth <= 1) iwidth = 2;
  
  /* -------------------------- */
  /* fit sliding triangle       */

  for (i = 0 ; i < n ; i++) {
  
    int j;
    float sum, fitval;

    sum = 0.;
    
    for (j = 0 ; j < iwidth ; j++) {
      int k = (j < iwidth/2 ? j : iwidth - j);
      float val = x[(i+j) % n] - (2. * (float) k / (float) iwidth * height);
      sum += val * val;
    }

    if ((i == 0) || (sum < fitval)) {
      fitval = sum;
      fitpos = i;
    }
  }  

  fitpos = (fitpos + (iwidth /2)) % n;

  return fitpos;
}

/* ----------------------------------------------- */
/* zero signal                                     */
/* ----------------------------------------------- */

static void sZero(float *x, int from, int to, int n) {
  int i;

  for (i = from ; i < to ; i++) {
    x[(i+n) % n] = 0.;
  }
}

/* ----------------------------------------------- */
/* integrate signal                                */
/* ----------------------------------------------- */

static float sIntegrate(float *x, int from, int to, int n) {
  int i;
  float sum;

  sum = 0.;
  for (i = from ; i < to ; i++) {
    sum += x[(i+n) % n];
  }
  
  return sum;
}


/* ----------------------------------------------- */
/* main						   */
/* ----------------------------------------------- */

main() {

  int i, npt, nread, seqlength, window, iwidth;
  float sum, tsum, best;
  float *x, *xx;
  char buffer[BUFSIZ];

  /* --------------------------------------------- */
  /* first line : npt, seqlength, window           */
  
  fgets(buffer, sizeof(buffer), stdin);
  
  if (   (sscanf(buffer, "%d%d%d", &npt, &seqlength, &window) != 3)
      || ((npt <= 0) || (seqlength <= 0) || (window <= 0))) {
    fprintf(stderr, "bad header (npt seqlength window) %s\n", buffer);
    exit(1);
  }

  /* --------------------------------------------- */
  /* allocate buffer                 	           */
  
  if (! (x = NEWN(float, npt))) {
    MemoryErreur("main", 2);
  }

  if (! (xx = NEWN(float, npt))) {
    MemoryErreur("main", 2);
  }

  /* --------------------------------------------- */
  /* other lines : data                            */
  
  nread = 0;
  while (fscanf(stdin, "%f", x + nread) == 1) {
    if (nread > npt) {
      (void) sprintf(buffer, "too many points (declared=%d)", npt);
      Erreur(buffer, 3);
    }
    nread++;
  }
  
  if (nread != npt) {
    (void) sprintf(buffer, "number of points does not match : declared=%d actual=%d", npt, nread);
    Erreur(buffer, 0);
  }

  /* --------------------------------------------- */
  /* print header 				   */

  printf("seqlength : %d\n", seqlength);
  printf("window    : %d\n", window);

  /* --------------------------------------------- */
  /* fit	 				   */

  iwidth = (int) ((float) window * (float) npt / (float) seqlength);

  sum = sIntegrate(x, 0, npt, npt);

  for (i = 1 ; i < npt ; i++)
    xx[i] = x[i];


  for (i = 0 ; i < NBEST ; i++) {
  
    tsum = sLookBest(xx, npt, iwidth);

    best = sFitTriangle(xx, npt, iwidth, tsum / (iwidth /2));
    
    tsum = sIntegrate(x, best - iwidth/2, best + iwidth/2, npt);
    
    best = sFitTriangle(xx, npt, iwidth, tsum / (iwidth /2));
    
    tsum = sIntegrate(x, best - iwidth/2, best + iwidth/2, npt);

    sZero(xx, best - iwidth/2, best + iwidth/2, npt);

    printf("best[%d] : %.0f %.0f %.2f\n", i, best, 1 + best * (float) seqlength / (float) npt,
                                          tsum * 100. / sum);
  }
  
  exit(0);
}
