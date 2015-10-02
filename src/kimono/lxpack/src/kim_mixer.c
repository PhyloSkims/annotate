/* ==================================================== */
/* @file: kim_mixer.c                                   */
/* count mixing functions                               */
/* @history:                                            */
/* @+   Apr. 97 <AV> first draft                        */
/* ==================================================== */

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "kimono.h"

/* ----------------------------------------------- */
/* basic mixer					   */
/*                                                 */
/* this mixer just pool the two strands on both    */
/* sides                                           */
/* ----------------------------------------------- */

void BasicMixer(int left[][2], int right[][2],
		int *countLeft, int *countRight,
		int nbins) {
  int i;
  
  for (i = 0 ; i <= nbins ; i++)
    countLeft[i] = countRight[i] = 0;

  for (i = 0 ; i <= nbins ; i++) {
    countLeft[i]  += left[i][0] + left[i][1];    /* leftDirect  + leftReverse  */
    countRight[i] += right[i][0] + right[i][1];  /* rightDirect + rightReverse */
  }
}

/* ----------------------------------------------- */
/* cross mixer					   */
/*                                                 */
/* this mixer cross the two strands on both        */
/* sides                                           */
/* ----------------------------------------------- */

void CrossMixer(int left[][2], int right[][2],
		int *countLeft, int *countRight,
		int nbins) {
  int i;
  
  for (i = 0 ; i <= nbins ; i++)
    countLeft[i] = countRight[i] = 0;

  for (i = 0 ; i <= nbins ; i++) {
    countLeft[i]  += left[i][0] + right[i][1];  /* leftDirect  + rightReverse */
    countRight[i] += right[i][0] + left[i][1];  /* rightDirect + leftReverse  */
  }
}

