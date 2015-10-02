/* ==================================================== */
/* @file: kim_alphabet.c                                */
/* alphabet handling                                    */
/* @history:                                            */
/* @+   Apr. 97 <AV> first draft                        */
/* ==================================================== */

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "kimono.h"

/* some static alphabets */

static char sDNA[]  = KIM_DNA;
static char sAA[]   = KIM_AA;

/* -------------------------------------------- */
/* API                                          */
/* -------------------------------------------- */

/* -------------------------------------------- */
/* free alphabet                                */
/* -------------------------------------------- */

Alphabet * FreeAlphabet(Alphabet *alpha) {
  if (alpha) FREE(alpha);
  return NULL;
}

/* -------------------------------------------- */
/* new alphabet					*/
/* -------------------------------------------- */

Alphabet *NewAlphabet(char *symbols) {
  Alphabet *alpha;
  
  if (alpha = NEW(Alphabet)) {
    (void) strcpy(alpha->symbol, symbols);
    alpha->length = strlen(alpha->symbol);
  }

  return alpha;
}

/* -------------------------------------------- */
/* copy alphabet                                */
/* -------------------------------------------- */

Alphabet *CopyAlphabet(Alphabet *alpha) {
  return NewAlphabet(alpha->symbol);
}

/* -------------------------------------------- */
/* standard alphabets                           */
/* -------------------------------------------- */

/* DNA */

Alphabet *NucleicAlphabet() {
  return NewAlphabet(sDNA);
}

/* AA  */

Alphabet *AminoAlphabet() {
  return NewAlphabet(sAA);
}

/* ----------------------------------------------- */
/* index of symbol in alphabet                     */
/* return : index of symbol in alphabet            */
/*          -1 if not found                        */
/* ----------------------------------------------- */

int AlphaIndex(Alphabet *alpha, char symbol) {

  char *pos = strchr(alpha->symbol, UPPER(symbol));
  
  return (pos ? (int) (pos - alpha->symbol) : -1);
}
