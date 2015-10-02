/* ---------------------------------------------------------------- */
/* Copyright (c) Atelier de BioInformatique                         */
/* @file: kim_genetic.c                                             */
/* @desc: bio. sequences utility functions                          */
/*                                                                  */
/* @history:                                                        */
/* @+       <Gloup> : Jan 96 : quick draft                          */
/* ---------------------------------------------------------------- */

#include <ctype.h>
#include <stdio.h>
#include <string.h>

#include "kimono.h"

#define GENETIC_CODE_INSTANCE
#include "Genetic.h"

static char sNuc[]     = DNA_ALPHA;
static char sAnuc[]    = C_DNA_ALPHA;
static char sGenNuc[]  = GEN_NUC_ALPHA;
static char sGenPro[]  = GEN_PRO_ALPHA;
static int  sNucNum[5] = {0, 1, 2, 3, 3};

/* ---------------------------------------------------- */
/* @Function :  int BaseComplement                      */
/* Purpose : return DNA/RNA-Iupac base complement       */
/* ---------------------------------------------------- */
int BaseComplement(int ch)
{
    char *c;

    if (c = strchr(sNuc, ch))
        return sAnuc[(int) (c - sNuc)];
    else
        return ch;
}

/* ---------------------------------------------------- */
/* @Function :  char * SeqComplement                    */
/* Purpose : return sequence complement                 */
/* ---------------------------------------------------- */
char *SeqComplement(char *str)
{
    char *s;

    for (s = str ; *s ; s++)
        *s = BaseComplement(*s);

    return str;
}

/* ---------------------------------------------------- */
/* @Static : int * sTranslateCodon                      */
/* Purpose : translate codon -> aa                      */
/* see  bio_codon_translate                             */
/* ---------------------------------------------------- */
static int sTranslateCodon(char *codon, int *code)
{
    int  i, base, hash;
    char *p;
    
    for (i = hash = 0 ; i < 3 ; i++) {
        if ((p = strchr(sGenNuc, *codon++)) != NULL) {
            base = ((int) (p - sGenNuc)) / 2;
            hash = (hash * 4) + sNucNum[base];
        }
        else {
            hash = 64;  /* bad letter in codon  */
            break;      /* or incomplete codon  */
        }
    }
    return (int) sGenPro[code[hash]];
}

/* ---------------------------------------------------- */
/* @Function :  int * CodonTranslate                    */
/* Purpose : return amino-acid                          */
/* input: codon char* 3 bases (in GEN_NUC_ALPHA)        */
/*        codid int   (see Genetic.h)                   */
/* output: aa in one letter code (in GEN_PRO_ALPHA)     */
/* ---------------------------------------------------- */
int CodonTranslate(char *codon, int codid)
{
    return sTranslateCodon(codon, theGeneticCode[codid].code);
}

/* ---------------------------------------------------- */
/* @Function :  int * SeqTranslate                      */
/* Purpose : translate sequence to protein              */
/* ---------------------------------------------------- */
char* SeqTranslate(char *seq, int codid)
{
    int  *code;
    char *ps, *ns;
    
    if ((codid < 0) || (codid >= GEN_MAX_CODES))
        return NULL;

    code = theGeneticCode[codid].code;

    for (ns = ps = seq ; ns[0] && ns[1] && ns[2] ; ns += 3)
        *ps++ = sTranslateCodon(ns, code);
    
    *ps = '\000';

    return seq;    
} 

