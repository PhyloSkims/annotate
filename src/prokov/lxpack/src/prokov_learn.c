/* ---------------------------------------------------------------- */
/* Copyright (c) Atelier de BioInformatique                         */
/* @file: prokov_learn.c                                            */
/* @desc: prokov - phase apprentissage                              */
/* @+           A partir de N sequences d'ADN (format fasta)        */
/* @+           calcule les occurences de tous les kupleCs sur les  */ 
/* @+           phases de lecture et ecrit les nombres d'occurence  */
/* @+           dans le fichier "matrice"                           */
/*                                                                  */
/* @history:                                                        */
/* @+   <Fred Nikitin + Marc Heuveline> : Aug 99 : first version    */
/* @+   <Gloup> : Oct 99 : last revised version                     */
/* @+   <Gloup> : Nov 99 : added Markov model for non-coding        */
/* @+   <Gloup> : Dec 00 : HelixWare port - V1.2                    */
/* ---------------------------------------------------------------- */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#ifdef HAS_GETOPT_H
#include HAS_GETOPT_H
#endif

#include "fasta_io.h"
#include "libkov.h"

#ifndef FILENAME_MAX
#define FILENAME_MAX 1024
#endif

#ifndef BUFSIZ
#define BUFSIZ 4096
#endif

/* -------------------------------------------- */
/* getopt globals                               */
/* -------------------------------------------- */

        extern char     *optarg;
        extern int      optind;

/* ---------------------------------------------------- */
/* copy src into dst                                    */
/* same as strcpy but works if src and dst overlaps     */
/* ---------------------------------------------------- */
static char *sStrCpy(char *dst, char *src)
{
    char *d;
    
    for (d = dst ; *src ; src++)
        *d++ = *src;

    *d = '\000';
    return dst;
}
    
/* ---------------------------------------------------- */
/* initialise a 0 un tableau d'occurences occ           */
/* ---------------------------------------------------- */
static void sInitOcc(int occ[POW_KMAX], int kuple)
{
    int j, powk;

    powk = Pow4(kuple);

    for (j = 0 ; j < powk ; j++)
        occ[j] = 0;
}

/* ---------------------------------------------------- */
/* calcule les occurences de tous les kuples dans seq   */
/* et met occ a jour                                    */
/*                                                      */
/* version 'en phase' pour codant                       */
/*                                                      */
/* ---------------------------------------------------- */
static void sAjoutOccC(char* seq, int occ[3][POW_KMAX], int kuple)
{
   int i, imax, phase, index;

   imax = strlen(seq) - kuple + 1;

   index = Codage(seq + 0, kuple);   /* valeur initiale */
   
   occ[0][index]++;     

   for (i = 1 ; i < imax ; i++) {
       index = SuiteCodage(seq + i, index, kuple);
       phase = i % 3;
       occ[phase][index]++;
   }   
}

/* ---------------------------------------------------- */
/* calcule les occurences de tous les kuples dans seq   */
/* et met occ a jour                                    */
/*                                                      */
/* version 'hors phase' pour non-codant                 */
/*                                                      */
/* ---------------------------------------------------- */
static void sAjoutOccN(char* seq, int occ[POW_KMAX], int kuple)
{
   int i, imax, index;

   imax = strlen(seq) - kuple + 1;

   index = Codage(seq + 0, kuple);   /* valeur initiale */
   
   occ[index]++;     

   for (i = 1 ; i < imax ; i++) {
       index = SuiteCodage(seq + i, index, kuple);
       occ[index]++;
   }   

}

/* ---------------------------------------------------- */
/* renormalise un tableau d'occurence                   */
/* utilise pour renomaliser le non-codant calcule sur   */
/* les 6 phases                                         */
/* ---------------------------------------------------- */
static void sNormalizeOcc(int occ[POW_KMAX], int kuple, int norm)
{
   int i, powk;
   
   powk = Pow4(kuple);
   
   for (i = 0 ; i < powk ; i++)
        occ[i] /= norm;
}

/* ---------------------------------------------------- */
/* occurences simulees pour le non-codant               */
/* somme du codant sur toutes les phases                */
/* == codant shuffle en mononucleotides                 */
/*                                                      */
/* ceci est utlise lorsqu'on ne connait pas d'exemples  */
/* de non-codant                                        */
/*                                                      */
/* ---------------------------------------------------- */
static void sNoCodingRan(int occpos[3][POW_KMAX], int occneg[3][POW_KMAX],
                         int occnoc[POW_KMAX], int kuple)
{
   int i, powk, phase, noc;

   powk = Pow4(kuple);

   for (i = 0 ; i < powk ; i++) {

       noc = 0;
       
       for (phase = 0 ; phase < 3 ; phase++) {
           noc += occpos[phase][i];
           noc += occneg[phase][i];
       }
       
       occnoc[i] = noc;
   }
}

/* ----------------------------------------------- */
/* printout help                                   */
/* ----------------------------------------------- */

#define PP (void) fprintf(stdout, 

static void sPrintHelp()
{
    PP  "------------------------------------------     \n");
    PP  " Prokov Version %s\n", VERSION);
    PP  "------------------------------------------     \n");
    PP  "synopsis :                                     \n");
    PP  "  Markov/Bayes based gene prediction           \n");
    PP  "  * learning phase *                           \n");
    PP  "                                               \n");
    PP  "usage: prokov_learn [options] [cod [nocod]]    \n");
    PP  "------------------------------------------     \n");
    PP  "options:                                       \n");
    PP  "                                               \n");      
    PP  "-a : write matrix in [A]scii format            \n");
    PP  "        default: write in binary format        \n");
    PP  "                                               \n");      
    PP  "-h : this [H]elp                               \n");
    PP  "                                               \n");      
    PP  "-i : [i]gnore first codon of each sequence     \n");
    PP  "        default: off                           \n");
    PP  "                                               \n");      
    PP  "-I : [I]gnore last codon of each sequence      \n");
    PP  "        default: off                           \n");
    PP  "                                               \n");      
    PP  "-k nn : use [K]uple of size nn                 \n");
    PP  "        for coding                             \n");
    PP  "        default: -k %d\n", KUP_DFT);
    PP  "                                               \n");      
    PP  "-K nn : use [K]uple of size nn                 \n");
    PP  "        for non-coding                         \n");
    PP  "        default: use mean of coding on 6 phases\n");
    PP  "                                               \n");      
    PP  "-o name : use name as [O]utput filename        \n");
    PP  "        default: -o %s\n", MATNAME);
    PP  "                                               \n");      
    PP  "-v : set [V]erbose mode                        \n");
    PP  "                                               \n");      
    PP  "------------------------------------------     \n");
    PP  " cod : file containing training genes          \n");
    PP  "             in Fasta format                   \n");
    PP  " note: genes may contain Start codon but not   \n");
    PP  "       Stop codon (see -I option)              \n");
    PP  " nocod : optional file containing noncoding    \n");
    PP  "             sequences in Fasta format         \n");
    PP  "------------------------------------------     \n");
}

#undef PP

/* ----------------------------------------------- */
/* printout usage and exit                         */
/* ----------------------------------------------- */

#define PP (void) fprintf(stderr, 

static int sExitUsage(int stat)
{
    PP  "usage: prokov_learn [options] [cod [nocod]]    \n");
    PP  "type \"prokov_learn -h\" for help              \n");

    if (stat)
        exit(stat);
            
    return stat;
}

#undef PP

/* ---------------------------------------------------- */
/*                PROGRAMME PRINCIPAL                   */
/* ---------------------------------------------------- */
int main(int argn, char *argv[]) 
{
    int  i, nbseqC, nbseqN;
    int  carg, errflag, aflag, vflag, giflag, iflag, rflag, kflag; 

    MarkovMatrix *pmat;

    FastaSequence *seq;

    FILE *filou;

    char buffer[BUFSIZ];
    char fileN[FILENAME_MAX], fileC[FILENAME_MAX],
         outname[FILENAME_MAX];

    /* ------------------------ */
    /* allocate matrix          */
    /* (avoid stack overflow    */
    /*  on some machines)       */
    /* ------------------------ */

    if (! (pmat = NEW(MarkovMatrix)))
        Erreur("Not enough memory", 10);

    /* ------------------------ */
    /* lecture des arguments    */
    /* ------------------------ */

    /* ------------------------ */
    /* valeurs par defaut       */
    
    errflag = 0;
 
    aflag  = 0;         /* binary output        */
    vflag  = 0;         /* not verbose          */
    iflag  = 0;         /* use first            */
    giflag = 0;         /* use last             */
    rflag  = 1;         /* use randonm coding   */
    kflag  = 0;         /* K unset              */

    pmat->kupleC = KUP_DFT;     /* kupleC = 3           */
    pmat->kupleN = KUP_DFT;     /* kupleN = 3           */

    (void) strcpy(pmat->version, VERSION);
    (void) strcpy(outname, MATNAME);

    while ((carg = getopt(argn, argv, "ahiIk:K:o:v")) != -1) {
            
        switch(carg) {

            case 'a' :              /* ascii mode   */
                aflag = 1;
                break;

            case 'h' :              /* help         */
                sPrintHelp();
                exit(0);
                break;

            case 'i' :              /* ignore first  */
                iflag = 1;
                break;

            case 'I' :              /* ignore last   */
                giflag = 1;
                break;

            case 'k' :              /* kupleC       */
                if (sscanf(optarg, "%d", &(pmat->kupleC)) != 1)
                   errflag++;
                break;

            case 'K' :              /* kupleN       */
                if (sscanf(optarg, "%d", &(pmat->kupleN)) != 1)
                   errflag++;
                rflag = 0;
                kflag = 1;
                break;

            case 'o' :              /* output name  */
                (void) strcpy(outname, optarg);
                break;

            case 'v' :              /* verbose      */
                vflag = 1;
                break;

            case '?' :              /* misusage     */
                errflag++;
        }
    }
    
    /* ------------------------ */
    /* accepts 0, 1 or 2        */
    /* arguments                */

    if ((argn - optind) == 0) {         /* 0 arg -> <stdin>     */
        rflag = 1;                      /* use random nocoding  */
        (void) strcpy(fileC, "<stdin>");
        (void) strcpy(fileN, "<none>");
    }
    else if ((argn - optind) == 1) {    /* 1 arg -> cod         */
        rflag = 1;                      /* use random nocoding  */
        (void) strcpy(fileC, argv[optind]);
        (void) strcpy(fileN, "<none>");
        if (! AssignToStdin(fileC))
           errflag++;
    }
    else if ((argn - optind) == 2) {    /* 2 args -> cod nocod  */
        rflag = 0;                      /* use real nocoding    */
        (void) strcpy(fileC, argv[optind]);
        (void) strcpy(fileN, argv[optind+1]);
        if (! AssignToStdin(fileC))
           errflag++;
        if (! AssignToStdin(fileN))
           errflag++;
    }
    else
        errflag++;

    /* ------------------------ */
    /* check arguments          */
    
    if ((rflag) || (! kflag))
        pmat->kupleN = pmat->kupleC;

    if (pmat->kupleC > K_MAX) {
        Erreur("kupleC is too large", 0);
        errflag++;
    }

    if (pmat->kupleN > K_MAX) {
        Erreur("kupleN is too large", 0);
        errflag++;
    }

    if (! (filou = OpenFile(outname, (aflag ? "w" : "wb"))))
        errflag++;

    /* ------------------------ */
    /* exit on usage error      */

    if (errflag) 
        (void) sExitUsage(1);

    /* ------------------------ */
    /* initialisations          */
    /* ------------------------ */

    pmat->powkC = Pow4(pmat->kupleC);
    pmat->powkN = Pow4(pmat->kupleN);

    for (i = 0 ; i < 3 ; i++) {
        sInitOcc(pmat->occpos[i], pmat->kupleC);
        sInitOcc(pmat->occneg[i], pmat->kupleC);
    }
    sInitOcc(pmat->occnc, pmat->kupleN);

    seq = NewFastaSequence();
    
    /* ------------------------ */
    /* calcul pour le  codant   */
    /* ------------------------ */

    nbseqC = 0;

    (void) AssignToStdin(fileC);

    while(ReadFastaSequence(stdin, seq)) {

        if (! seq->ok) {
           (void) sprintf(buffer, "invalid sequence %s", seq->name);
           Erreur(buffer, 0);
           continue;
        }
         
        nbseqC++;
        
        if (vflag)
           (void) fprintf(stderr, "Gene: %s [%ld b]\n", seq->name, seq->length);
           
        (void) RNAtoDNA(Upper(seq->seq));

        if (iflag && (seq->length >= 3)) {      /* remove start */
           (void) sStrCpy(seq->seq, seq->seq+3);
           seq->length -= 3;
        }
        
        if (giflag && (seq->length > 3)) {      /* remove stop  */
           seq->seq[seq->length-4] = '\000';
           seq->length -= 3;
        }
        
        sAjoutOccC(seq->seq, pmat->occpos, pmat->kupleC);       /* brin direct          */
                                                                /* seulement            */
#if 0   
        Complement(Reverse(seq->seq));
        sAjoutOccC(seq->seq, pmat->occneg, pmat->kupleC);       /* brin complementaire  */
                                                                /* not needed in V1.2   */
#endif

    }

                                                                /* brin complementaire  */
    ComputeNegOccurences(pmat->occneg, pmat->occpos, pmat->kupleC);

    /* --------------------------- */
    /* calcul pour le noncodant    */
    /* --------------------------- */
   
    nbseqN = 0;

    if (rflag) {
       sNoCodingRan(pmat->occpos, pmat->occneg, pmat->occnc, pmat->kupleC);
    }
    else {
   
       (void) AssignToStdin(fileN);
   
       while(ReadFastaSequence(stdin, seq)) {

           if (! seq->ok)
              continue;

           nbseqN++;

           if (vflag)
              (void) fprintf(stderr, "Non-Coding: %s [%ld b]\n",
                                     seq->name, seq->length);

           (void) RNAtoDNA(Upper(seq->seq));

           sAjoutOccN(seq->seq, pmat->occnc, pmat->kupleN);  /* brin direct     */

           Complement(Reverse(seq->seq));

           sAjoutOccN(seq->seq, pmat->occnc, pmat->kupleN);  /* complementaire  */
       }
    }
   
    sNormalizeOcc(pmat->occnc, pmat->kupleN, 6);
    
    /* --------------------------- */
    /* ecriture matrice            */
    /* --------------------------- */

    if (! WriteMatrix(pmat, aflag, filou))
       exit(3);

    (void) fclose(filou);

    /* --------------------------- */
    /* fin de programme            */
    /* --------------------------- */

    FREE(pmat);

    FreeFastaSequence(seq);

    return 0 ;
}
