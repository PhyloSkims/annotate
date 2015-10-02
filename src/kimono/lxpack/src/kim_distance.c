/* ==================================================== */
/* @file: kim_distance.c                                */
/* general utilities functions                          */
/* @history:                                            */
/* @+   Apr. 97 <AV> first draft                        */
/* ==================================================== */

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "kimono.h"

#define SQRTF(x) (float) sqrt((float) x)
#define FABSF(x) (float) fabs((float) x)

/* ---------------------- */
/* API                    */
/* ---------------------- */

/*
 * for all distance functions :
 * int bins1[], int bins2[] : count distributions
 * int nbins                : total number of bins
 * int ibin		   : if < 0  : compute distance using all bins
 *                         : if >= 0 : compute distance only for ibin (partial distance) 
 */

/* ---------------------- */
/* chisquare distance     */

float ChisTwo(int bins1[], int bins2[], int nbins, int ibin)
{
  int j, jmin, jmax;
  float num, s1, s2, c1, c2, sj, dist;

  s1 = s2 = 0.;
  for (j = 0 ; j < nbins ; j++) {
    s1 += bins1[j];
    s2 += bins2[j];
  }
  
  if (s1 == 0.) s1 = 1.0;
  if (s2 == 0.) s2 = 1.0;
  
  c1 = SQRTF(s2/s1);
  c2 = SQRTF(s1/s2);

  jmin = (ibin < 0 ? 0 : ibin);
  jmax = (ibin < 0 ? nbins : ibin + 1);

  dist = 0.0;
  
  for (j = jmin ; j < jmax ; j++) {
    sj = bins1[j] + bins2[j];
    if (sj == 0.) sj = 1.;
    num = (c1 * (float) bins1[j]) - (c2 * (float) bins2[j]);
    dist += (num * num) / sj;
  }
  
  return SQRTF(dist);
}

/* ---------------------- */
/* euclidian distance     */

float Euclid(int bins1[], int bins2[], int nbins, int ibin)
{
  int j, jmin, jmax;
  float num, s1, s2, dist;

  s1 = s2 = 0.;
  for (j = 0 ; j < nbins ; j++) {
    s1 += bins1[j];
    s2 += bins2[j];
  }
  
  if (s1 == 0.) s1 = 1.0;
  if (s2 == 0.) s2 = 1.0;

  jmin = (ibin < 0 ? 0 : ibin);
  jmax = (ibin < 0 ? nbins : ibin + 1);

  dist = 0.0;

  for (j = jmin ; j < jmax ; j++) {
    num = ((float) bins1[j]/s1) - ((float) bins2[j]/s2);
    dist += num*num;
  }
  
  return SQRTF(dist);
}

/* ---------------------- */
/* hellinger distance     */

float HellTwo(int bins1[], int bins2[], int nbins, int ibin)
{
  int j, jmin, jmax;
  float num, s1, s2, dist;

  s1 = s2 = 0.;
  for (j = 0 ; j < nbins ; j++) {
    s1 += bins1[j];
    s2 += bins2[j];
  }
  
  if (s1 == 0.) s1 = 1.0;
  if (s2 == 0.) s2 = 1.0;

  jmin = (ibin < 0 ? 0 : ibin);
  jmax = (ibin < 0 ? nbins : ibin + 1);

  dist = 0.0;

  for (j = jmin ; j < jmax ; j++) {
    num = SQRTF((float) bins1[j]/s1) - SQRTF((float) bins2[j]/s2);
    dist += num*num;
  }
  
  return SQRTF(dist);
}

/* ----------------------------------------------- */
/* compute skew sign for two-binned data           */
/* ----------------------------------------------- */

int SignTwoBins(int bins1[], int bins2[]) {

  int s1, s2;
  
  s1 = bins1[0] + bins1[1];
  s2 = bins2[0] + bins2[1];
  
  if (s1 == 0) s1 = 1;
  if (s2 == 0) s2 = 1;

  return ( (((float) bins1[0] / (float) s1) >= ((float) bins2[0] / (float) s2))
           ? 1 : -1);
}
