/* ==================================================== */
/* @file: kim_indexer.c                                 */
/* indexing functions                                   */
/* @history:                                            */
/* @+   Jul. 06 <AV> first draft                        */
/* ==================================================== */

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "kimono.h"

/* -------------------------------------------- */
/* internal functions                           */
/* -------------------------------------------- */

/* ----------------------------------------------- */
/* codon index                                     */
/* codon: start of codon sequence                  */
/*        read left to right for direct strand     */
/*             right to left for reverse strand    */
/* ----------------------------------------------- */

static int sCodonIndex(char *codon, Alphabet *alpha, int strand)
{
  int i, code;

  code = 0;
    
  for (i = 0 ; i < 3 ; i++) {    /* loop on codon positions */
  
    char symb;
    
    if (strand > 0) {
      symb = codon[i];                   /* move to left    */
    }
    else {
      symb = codon[-i];                  /* move to right   */
      symb = BaseComplement(symb);       /* complement      */
    }

    int k = AlphaIndex(alpha, symb);
    
    if (k < 0) 
      return -1;                         /* invalid symbol  */
    
    code = (code << SHIFT_NUC) | k;      /* encode codon    */
  }

  return code;  
}

/* ----------------------------------------------- */
/* basic sequence indexing                         */
/* ----------------------------------------------- */

static void sBasicIndexer(Sequence *seq) {
  int i;
  for (i = 0 ; i < seq->fasta->length ; i++) {
    seq->cell[i].index = AlphaIndex(seq->alphabet, seq->fasta->seq[i]);
  }
}

/* ----------------------------------------------- */
/* reset mask                                      */
/* ----------------------------------------------- */

static void sSetMask(Sequence *seq, int mask) {
  int i;
  for (i = 0 ; i < seq->fasta->length ; i++) {
    seq->cell[i].mask = mask;
  }
}

/* ----------------------------------------------- */
/* mask masked symbol positions                    */
/* ----------------------------------------------- */

static void sMaskMaskedSymbol(Sequence *seq) {
  int i;
  for (i = 0 ; i < seq->fasta->length ; i++) {
    if (seq->fasta->seq[i] == MASKED_SYMBOL)
      seq->cell[i].mask |= EXCLUDE_MASK;
  }
}

/* ----------------------------------------------- */
/* mask overlaps                                   */
/* there is an overlap when a position is both on  */
/* direct and reverse strand.                      */
/* this situation makes trouble since the index is */
/* wrong on one strand.                            */
/* note [AV] : overlaps on the same strand are not */
/*             detected yet but this situation is  */
/*             not such a big trouble.             */
/* ----------------------------------------------- */

static void sMaskOverlap(Sequence *seq) {
  int i;
  for (i = 0 ; i < seq->fasta->length ; i++) {
    if (seq->cell[i].mask & (DIRECT_MASK & REVERSE_MASK))
      seq->cell[i].mask |= EXCLUDE_MASK;
  }
}


/* -------------------------------------------- */
/* API                                          */
/* -------------------------------------------- */

/* ----------------------------------------------- */
/* no zone indexing                                */
/* (no strand, no phase)                           */
/* ----------------------------------------------- */

void NoZoneIndexer(Sequence *seq, SetOfIntervals *zones, int phase) {
/*ARGSUSED*/

  /* ---------------------------------- */
  /* assume everything is direct        */
  
  sSetMask(seq, DIRECT_MASK);

  /* ---------------------------------- */
  /* exclude masked symbol              */
  
  sMaskMaskedSymbol(seq);

  /* ---------------------------------- */
  /* just use basic indexer             */
  
  sBasicIndexer(seq);
}

/* ----------------------------------------------- */
/* zone exclusion indexing                         */
/* (no strand, no phase)                           */
/* ----------------------------------------------- */

void ExcludeIndexer(Sequence *seq, SetOfIntervals *zones, int phase) {
/*ARGSUSED*/

  int      i, j;
  Interval inter;

  /* ---------------------------------- */
  /* assume everything is direct        */
  
  sSetMask(seq, DIRECT_MASK);

  /* ---------------------------------- */
  /* mask zone positions                */

  for (i = 0; i < zones->size ; i++) {
    inter = zones->interval[i];
    for (j = inter.from ; j <= inter.to ; j++)
      seq->cell[j].mask |= EXCLUDE_MASK;
  }

  /* ---------------------------------- */
  /* exclude masked symbol              */
  
  sMaskMaskedSymbol(seq);
  
  /* ---------------------------------- */
  /* then use basic indexer             */

  sBasicIndexer(seq);
}

/* ----------------------------------------------- */
/* gene indexing                                   */
/* (strand, no phase)                              */
/* note: alpha is unused                           */
/* ----------------------------------------------- */

void GeneIndexer(Sequence *seq, SetOfIntervals *zones, int phase) {
/*ARGSUSED*/

  int i;

  /* ---------------------------------- */
  /* first reset mask, excluding all    */

  sSetMask(seq, EXCLUDE_MASK);

  /* ---------------------------------- */
  /* then unmask start of genes         */
  
  for (i = 0; i < zones->size ; i++) {
    Interval inter = zones->interval[i];
    int pos  = (inter.strand > 0 ? inter.from : inter.to);
    int mask = (inter.strand > 0 ? DIRECT_MASK : REVERSE_MASK);
    seq->cell[pos].mask = (seq->cell[pos].mask & INCLUDE_MASK) | mask;
    seq->cell[pos].index = (inter.strand > 0 ? 0 : 1);
  }

  /* ---------------------------------- */
  /* finally, exclude overlaps          */
  /* when two genes start at same pos   */

  sMaskOverlap(seq);
}

/* ----------------------------------------------- */
/* phased nucleotide indexing                      */
/* phase: binary mask indicating valid phases      */
/*        1[01] : 1 ; 2[10] : 2 ; 4[100] : 3       */
/*        or any combination eg: 3[11] : 12        */
/*                               7[111] : 123      */
/* ----------------------------------------------- */
void PhaseIndexer(Sequence *seq, SetOfIntervals *zones, int phase) {
/*ARGSUSED*/

  int i, j;

  /* ---------------------------------- */
  /* first reset mask, excluding all    */

  sSetMask(seq, EXCLUDE_MASK);

  /* ---------------------------------- */
  /* then unmask valid positions        */
  
  for (i = 0; i < zones->size ; i++) {
    Interval inter = zones->interval[i];
    int iphase = 0;
    if (inter.strand > 0) {
      for (j = inter.from ; j <= inter.to ; j++) {
        if (((0x1 << iphase) & phase) != 0) 
          seq->cell[j].mask = (seq->cell[j].mask & INCLUDE_MASK) | DIRECT_MASK;
        iphase = (iphase + 1) % 3;
      }
    }
    else {
      for (j = inter.to ; j >= inter.from ; j--) {
        if (((0x1 << iphase) & phase) != 0)
          seq->cell[j].mask = (seq->cell[j].mask & INCLUDE_MASK) | REVERSE_MASK;
        iphase = (iphase + 1) % 3;
      }
    }
  }
  
  /* ---------------------------------- */
  /* exclude masked symbol              */
  
  sMaskMaskedSymbol(seq);

  /* ---------------------------------- */
  /* exclude overlaps                   */

  sMaskOverlap(seq);
  
  /* ------------------------------------------ */
  /* finally, use basic indexer                 */

  sBasicIndexer(seq);
}

/* ----------------------------------------------- */
/* codon indexing                                  */
/* ----------------------------------------------- */

void CodonIndexer(Sequence *seq, SetOfIntervals *zones, int phase) {
  int i, j;

  /* ---------------------------------- */
  /* reset mask, excluding all          */

  sSetMask(seq, EXCLUDE_MASK);

  /* ------------------------------------------ */
  /* now unmask phase 1 and compute codon index */

  for (i = 0; i < zones->size ; i++) {
    Interval inter = zones->interval[i];
    int iphase = 0;
    if (inter.strand > 0) {
      for (j = inter.from ; j <= inter.to - 2 ; j++) { // ZZ to - 2
        if (iphase == 0) {
          int k = sCodonIndex(seq->fasta->seq + j, seq->alphabet, 1);
          seq->cell[j].mask = (seq->cell[j].mask & INCLUDE_MASK) | DIRECT_MASK;
          seq->cell[j].index = k;
        }
        iphase = (iphase + 1) % 3;
      }
    }
    else {
      for (j = inter.to ; j >= inter.from + 2 ; j--) { // ZZ from + 2
        if (iphase == 0) {
          int k = sCodonIndex(seq->fasta->seq + j, seq->alphabet, -1);
          seq->cell[j].mask = (seq->cell[j].mask & INCLUDE_MASK) | REVERSE_MASK;
          seq->cell[j].index = k;
        }
        iphase = (iphase + 1) % 3;
      }
    }
  }

  /* ---------------------------------- */
  /* exclude overlaps                   */

  sMaskOverlap(seq);

}


