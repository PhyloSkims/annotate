/* ==================================================== */
/* @file: kim_counter.c                                 */
/* counting functions                                   */
/* @history:                                            */
/* @+   Jul. 06 <AV> first draft                        */
/* ==================================================== */

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "Genetic.h"
#include "kimono.h"

/* holds nucleotide enrichment for each codon */

static float sCodonScore[NB_NUC][NB_CODONS];
static int sScoreInited = 0;

/* hold DNA amhabet

static char sDNA[]  = KIM_DNA;


/* -------------------------------------------- */
/* API                                          */
/* -------------------------------------------- */

/* ----------------------------------------------- */
/* basic nucleotide counting                       */
/*                                                 */
/* count : second dimension is for strand          */
/* increment : 1, -1 or 0 (0 means reset then 1)   */
/* ----------------------------------------------- */

void BasicCounter(Sequence *seq, int from, int to, 
                  int count[][2], int nbins, int increment) {

  int i, len;
  
  if (increment == 0) {        /* this is a reset */
    for (i = 0 ; i <= nbins ; i++)
      count[i][0] = count[i][1] = 0;
    increment = 1;
  }
  
  len = seq->fasta->length;

  for (i = from ; i < to ; i++) {

    int mask = seq->cell[i % len].mask;
    if (mask & EXCLUDE_MASK)
      continue;

    int index = seq->cell[i % len].index;
    if ((index < 0) || (index > nbins)) /* other symbol */
      index = nbins;

    if (mask & DIRECT_MASK)
      count[index][0] += increment;
    if (mask & REVERSE_MASK)
      count[index][1] += increment;
  }
}

/* ----------------------------------------------- */
/* counting skewed codon                           */
/* note: don't forget to init counter before       */
/* InitCodonSkewCounter(code)                      */
/* ----------------------------------------------- */

void InitCodonSkewCounter(int code) {
    CodonSkew(code, sCodonScore);
    sScoreInited = 1;
}

void CodonSkewCounter(Sequence *seq, int from, int to,
                      int count[][2], int nbins, int increment) {

  int i, len;

  if (! sScoreInited) {
    Erreur("CodonSkewCounter: initialized with default genetic code", 0); 
    InitCodonSkewCounter(GEN_CODE_UNIVL);
  }
  
  if (increment == 0) {        /* this is a reset */
    for (i = 0 ; i <= nbins ; i++)
      count[i][0] = count[i][1] = 0;
    increment = 1;
  }
  
  len = seq->fasta->length;

  for (i = from ; i < to ; i++) {
  
    int mask = seq->cell[i % len].mask;
    if (mask & EXCLUDE_MASK)
      continue;

    int index = seq->cell[i % len].index;    /* codon index   */
    if ((index < 0) || (index > NB_CODONS))  /* invalid codon */
      continue;

    int strand = ((mask & DIRECT_MASK) ? 0 : 1);

    int nuc;
    for (nuc = 0 ; nuc < NB_NUC ; nuc++) {

      float score = sCodonScore[nuc][index]; /* enrichment for symbol nuc */

      /* ---------------------------------------------------------------- */
      /* important property :						  */
      /* if codon causes enrichment/depletion in X and if the gene is     */
      /* located on the reverse strand then choice of this codon causes   */
      /* enrichment/depletion in ~X on the direct strand.		  */
      /* Therefore for reverse strand codons we should increment the bin  */
      /* associated to ~X						  */
      /*                                                                  */
      
      int bin = nuc;
      
      if (mask & REVERSE_MASK) {
        char compl = BaseComplement(seq->alphabet->symbol[nuc]);
        bin = AlphaIndex(seq->alphabet, compl);
      }

      /* ---------------------------------------------------------------- */
      /* increment bin associated to null, enrichment, depletion	  */

      bin = 3 * bin + ((score == 0.) ? 0 : ((score > 0) ? 1 : 2));

      count[bin][strand] += increment;
    }

    count[3 * NB_NUC][strand] += increment;             /* total # codons */
  }

}
