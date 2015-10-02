/* ==================================================== */
/* @file: kim_sequence.c                                */
/* sequence handling                                    */
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
static int sReallocCells(Sequence *seq) {
  SeqCell *buf;
  
  buf = (   seq->cell
          ? REALLOC(SeqCell, seq->cell, seq->fasta->length)
          : NEWN(SeqCell, seq->fasta->length));

  if (! buf) {
    MemoryErreur("NewSequence", 2);
    return 0;
  }
  
  seq->cell = buf;
  return 1;
}

/* -------------------------------------------- */
/* API                                          */
/* -------------------------------------------- */

/* -------------------------------------------- */
/* Sequence deallocation                        */
/* -------------------------------------------- */

Sequence *FreeSequence(Sequence *seq) {
  if (seq) {
    if (seq->alphabet) FreeAlphabet(seq->alphabet);
    if (seq->fasta)    FreeFastaSequence(seq->fasta);
    if (seq->cell)     FREE(seq->cell);
    FREE(seq);
  }
  return NULL;
}

/* -------------------------------------------- */
/* Sequence allocation                          */
/* -------------------------------------------- */

Sequence *NewSequence(Alphabet *alpha) {
  Sequence *seq;
  
  if (! (seq = NEW(Sequence)))
    MemoryErreur("NewSequence", 2);
    
  seq->alphabet = NULL;
  seq->fasta    = NULL;
  seq->cell     = NULL;
  
  if (! (seq->alphabet = CopyAlphabet(alpha)))
    MemoryErreur("NewSequence", 2);
  
  if (! (seq->fasta = NewFastaSequence()))
    MemoryErreur("NewSequence", 2);
  
  return seq;
}

/* -------------------------------------------- */
/* read sequence                                */
/* -------------------------------------------- */

int ReadSequence(FILE *streamin, Sequence *seq) {
  int ok;
  
  ok =    ReadFastaSequence(streamin, seq->fasta)
       && seq->fasta->ok
       && sReallocCells(seq);
       
  return ok;
}
