/* ---------------------------------------------------------------- */
/* Copyright (c) Atelier de BioInformatique                         */
/* @file: libkov_util.c                                             */
/* @desc: utilitaires prokov                                        */
/*                                                                  */
/* @history:                                                        */
/* @+   <Fred Nikitin + Marc Heuveline> : Aug 99 : first version    */
/* @+   <Gloup> : Oct 99 : last revised version                     */
/* @+   <Gloup> : Dec 00 : HelixWare port                           */
/* ---------------------------------------------------------------- */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "fasta_io.h"
#include "libkov.h"

#define BIO_FRAME       1       /* set to 1 for biological frames */

/* -------------------------------------------- */
/* error reporting                              */
/* -------------------------------------------- */
int Erreur(char *msg, int stat)
{
    (void) fprintf(stderr, "*Error* [%d] %s\n", stat, msg);

    (void) fflush(stderr);
        
    if (stat > 0)
       exit(stat);

    return stat;
}

/* ---------------------------------------------------- */
/* complemente la chaine str                            */
/* ---------------------------------------------------- */
char *Complement(char *str)
{
    char *s;
    
    for (s = str ; *s ; s++) {

        switch(*s) {
            case 'A' : *s = 'T'; break;  
            case 'T' : *s = 'A'; break;
            case 'C' : *s = 'G'; break;
            case 'G' : *s = 'C'; break;
        }

     }

     return str;   
}

/* ---------------------------------------------------- */
/* renverse la chaine str                               */
/* ---------------------------------------------------- */
char *Reverse(char *str)
{
    char *sb, *se, c;

    if (! str)
        return str;
            
    sb = str;
    se = str + strlen(str) - 1;

    while(sb <= se) {
        c    = *sb;
        *sb++ = *se;
        *se-- = c;
    }

    return str;
}

/* ---------------------------------------------------- */
/* uppercase la chaine str                              */
/* ---------------------------------------------------- */
char *Upper(char *str)
{
    char *s;
    
    for (s = str ; *s ; s++) 
        if (islower(*s))
           *s = toupper(*s);
           
    return str;
}

/* ---------------------------------------------------- */
/* convert RNA to DNA                                   */
/* ---------------------------------------------------- */
char *RNAtoDNA(char *str)
{
    char *s;
    
    for (s = str ; *s ; s++) 
        if (*s == 'U')
           *s = 'T';
           
    return str;
}

/* ---------------------------------------------------- */
/*  valeur des symboles                                 */
/* A=0; C=1; G=2; T=3 ; Other = 0                       */
/* ---------------------------------------------------- */
int Valeur(char s)
{
    switch (s) {
        case 'A' : return 0;
        case 'C' : return 1;
        case 'G' : return 2;
        case 'T' : return 3;
    }

    return 0; /* ? use A if unknown */
}

/* ---------------------------------------------------- */
/*  valeur du symbol complŽmentaire                     */
/* ---------------------------------------------------- */
int ValeurComplement (int valeur)
{
    switch (valeur) {
        case 0 : return 3; /* A -> T */ 
        case 1 : return 2; /* C -> G */
        case 2 : return 1; /* G -> C */
        case 3 : return 0; /* T -> A */
    }

     return 0; /* ? use A if unknown */ 
}

/* ---------------------------------------------------- */
/* calcule 4^n                                          */
/* tabulated version                                    */
/* ---------------------------------------------------- */
int Pow4(int n)
{
    int i, pow;

    static int sInited = 0;
    static int sTable[K_MAX];

    if (sInited)
       return sTable[n];
       
    pow = 1;

    for (i = 0 ; i < K_MAX ; i++) {
        sTable[i] = pow;
        pow *= 4;
    }
    
    sInited = 1;
    
    return sTable[n];
}

/* ---------------------------------------------------- */
/* hash coding d'un kuple (pour 4 symboles)             */
/* ---------------------------------------------------- */
int Codage(char *seq, int kuple)
{
    int i, index;

    index = 0;
    for (i = 0 ; i < kuple ; i++)
        index = (index << 2) | Valeur(seq[i]);

    return index;
}

/* ---------------------------------------------------- */
/* hash coding du complement-inverse du k-uple de       */
/* hashcode index                                       */
/* ---------------------------------------------------- */
int CodageComplement(int index, int kuple)
{
    int i, cind;

    cind = 0;
    for (i = 0 ; i < kuple ; i++) {
        cind = (cind << 2) |  ValeurComplement(index & 0x03);
        index = (index >> 2);
    }

    return cind;
}

/* ---------------------------------------------------- */
/* hash coding d'un kuple (pour 4 symboles) lorsque le  */
/* kuple precedent a deja ete calcule                   */
/* ---------------------------------------------------- */
int SuiteCodage(char *seq, int prec, int kuple)
{
    int new, mask;

    mask = Pow4(kuple) - 1; /* 0x7fffffff >> (32 - (2*k)) */

    new = ((prec << 2) | Valeur(seq[kuple-1])) & mask;

    return new;
}

/* ---------------------------------------------------- */
/* recherche du codon suivant dans une liste            */
/* ---------------------------------------------------- */
int FindCodon(char *seq, int lenseq, char *codons, int from, int to)
{
    int len;
    char save;
    char *sf, *st, *ss, *found;

    len = lenseq - 1;
    
    from = MIN(MAX(0, from), len);
    to   = MIN(MAX(from, to), len);

    /* round to multiple of 3 */

    len = ((to - from + 1) / 3) * 3; 
    to  = from + len - 1;
         
    sf = seq + from;
    st = seq + to;

    for (ss = sf ; ss <= st ; ss += 3) {
        save  = ss[3];
        ss[3] = '\000';
        found = strstr(codons, ss);
        ss[3] = save;
        if (found) 
           return ((int) (ss - sf)) + from;
    }

    return -1;
}

/* ---------------------------------------------------- */
/* calcul du frame de calcul correspondant au frame de  */
/* display dframe  - convention bio -                   */
/* pour la 'convention genemark' mettre BIO_FRAME a 0   */
/* note: les frames sont numerotes 0, 1, 2              */
/* ---------------------------------------------------- */
int InternalFrame(int dframe, int direct, int seqlen)
{
    int iframe;
    
    if (direct) {                               /* ---- direct strand  ---- */
       iframe = (3 - dframe) % 3;                       /* math convertion  */
    }
    else {                                      /* ---- reverse strand ---- */

#if BIO_FRAME
       iframe = (5 - dframe + ((seqlen+1) % 3)) % 3;    /* bio convention   */
#else
       iframe = dframe;                                 /* GM convention    */
#endif
       iframe = (3 - iframe) % 3;                       /* math convertion  */
    }

    return iframe;
}

/* ---------------------------------------------------- */
/* calcul des occurences en phase negative a partir     */
/* des occurences en phase positive                     */
/* ---------------------------------------------------- */

void ComputeNegOccurences(int occneg[3][POW_KMAX],
                          int occpos[3][POW_KMAX], int kuple)
{
    int  i, j, powk, phase, phaseinv;

    powk = Pow4(kuple);
   
    for (phase = 0 ; phase < 3 ; phase ++) {
        phaseinv = (3 - (kuple + phase) % 3 ) % 3; 
        for (i = 0 ; i < powk ; i++) {
            j = CodageComplement(i, kuple);
            occneg[phase][i] = occpos[phaseinv][j];
        }
    }
}

