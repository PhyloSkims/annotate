/* ==================================================== */
/* @file: kim_codonskew.c                               */
/* general utilities functions                          */
/* @history:                                            */
/* @+   Apr. 97 <AV> first draft                        */
/* ==================================================== */

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "kimono.h"

/* ----------------------------------------------- */
/* make string buffer of all codons :              */
/* "AAA/AAC/AAT/AAG/ACA/....../GGG/"               */
/* ----------------------------------------------- */

static void sMakeStdCodons(char *buffer, char *alphaN)
{
  int i, j, k, n;

  n = strlen(alphaN);
    
  for (i = 0 ; i < n ; i++)
  for (j = 0 ; j < n ; j++)
  for (k = 0 ; k < n ; k++) {
    *buffer++ = alphaN[i];
    *buffer++ = alphaN[j];
    *buffer++ = alphaN[k];
    *buffer++ = '/';
  }
  *buffer = '\000';  
}

/* ----------------------------------------------- */
/* make string buffer of synonymous codons coding  */
/* for aa                                          */
/* ----------------------------------------------- */

static void sMakeAaCodons(char *buffer, char *alphaN, int aa, int code)
{
  int     aai;
  char    *c, codons[4*NB_CODONS+1];
  
  sMakeStdCodons(codons, alphaN);  
  
  for (c = codons ; *c ; c += 4) {
    aai = CodonTranslate(c, code);
    if (aai == aa) {
      strncpy(buffer, c, 3);
      buffer[3] = '/';
      buffer += 4;
    }
  }
  *buffer = '\000';
}

/* ----------------------------------------------- */
/* return index of aa encoded by codon             */
/* or -1 if not found                              */
/* ----------------------------------------------- */

static int sAaIndex(char *codon, char *alphaA, int code)
{
  int  aa;
  char *paa;
  
  aa = CodonTranslate(codon, code);
  
  if ((paa = strchr(alphaA, aa)) != 0) 
    return (int) (paa - alphaA);
    
  return -1;
}

/* ----------------------------------------------- */
/* count # occurences of char mark in buffer       */
/* ----------------------------------------------- */

static int sCount(char *buffer, int mark)
{
  int count;
  
  for (count = 0 ; *buffer ; buffer++) {
    if (*buffer == mark)
      count++;
  }
  
  return count;
}

/* ----------------------------------------------- */
/* skew score of one codon                         */
/* ----------------------------------------------- */

static float sCodonSkew(char *codon, char *alphaN, char *alphaA, int symb, int code) {

  int iaa, nsyn, n0, nM;
  char synon[4*NB_CODONS+1], buf[4];
  
  strncpy(buf, codon, 3);
  buf[3] = '\000';

  iaa = sAaIndex(buf, alphaA, code);

  sMakeAaCodons(synon, alphaN, alphaA[iaa], code);
  
  n0 = sCount(buf, symb);  
  
  nsyn = strlen(synon) / 4;

  nM = sCount(synon, symb);

  return ((float) (n0) - ((float) (nM) / (float) nsyn));
}

/* ----------------------------------------------- */
/* API                                             */
/* ----------------------------------------------- */

/* ----------------------------------------------- */
/* skew scores of all codons                       */
/* ----------------------------------------------- */

void CodonSkew(int code, float score[NB_NUC][NB_CODONS]) {

    int  i, j;
    char *trip;
    char alpha[NB_NUC+1], codons[4*NB_CODONS+1];

    (void) strcpy(alpha, KIM_DNA);    

    sMakeStdCodons(codons, alpha);
    
    for (i = 0 ; i < NB_NUC ; i++) {
      for (j = 0, trip = codons ; j < NB_CODONS ; j++, trip += 4) {
        score[i][j] = sCodonSkew(trip, alpha, KIM_AA, alpha[i], code);
      }
    }
}
