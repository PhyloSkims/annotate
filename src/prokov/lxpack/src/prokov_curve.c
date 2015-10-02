/* ---------------------------------------------------------------- */
/* Copyright (c) Atelier de BioInformatique                         */
/* @file: prokov_curve.c                                            */
/* @desc: prokov - phase detection : trace courbes                  */
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

/* -------------------------------------------- */
/* getopt globals                               */
/* -------------------------------------------- */

        extern char     *optarg;
        extern int      optind;

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
    PP  "  * detection phase - curve *                  \n");
    PP  "                                               \n");
    PP  "usage: prokov_curve [options] seqfile          \n");
    PP  "------------------------------------------     \n");
    PP  "options:                                       \n");
    PP  "                                               \n");      
    PP  "-a : read matrix in [A]scii format             \n");
    PP  "        default: read in binary format         \n");
    PP  "                                               \n");      
    PP  "-h : this [H]elp                               \n");
    PP  "                                               \n");      
    PP  "-m name : [M]atrix filename                    \n");
    PP  "        default: -m %s\n", MATNAME);
    PP  "                                               \n");      
    PP  "-p nn : use nn %% as prior coding probability  \n");
    PP  "        default: -p %d\n", PCODANT);
    PP  "                                               \n");      
    PP  "-s nn : window [S]tep                          \n");
    PP  "        default: -s %d\n", WIN_STEP);
    PP  "                                               \n");      
    PP  "-v : set [V]erbose mode                        \n");
    PP  "                                               \n");      
    PP  "-w nn : window [W]idth                         \n");
    PP  "        default: -w %d\n", WIN_WIDTH);
    PP  "                                               \n");      
    PP  "------------------------------------------     \n");
    PP  " seqfile : file containing one sequence        \n");
    PP  "             in Fasta format                   \n");
    PP  "------------------------------------------     \n");
}

#undef PP

/* ----------------------------------------------- */
/* printout usage and exit                         */
/* ----------------------------------------------- */

#define PP (void) fprintf(stderr, 

static int sExitUsage(int stat)
{
    PP  "usage: prokov_curve [options] seqfile  \n");
    PP  "type \"prokov_curve -h\" for help      \n");

    if (stat)
        exit(stat);
            
    return stat;
}

#undef PP

/* ---------------------------------------------------- */
/*                PROGRAMME PRINCIPAL                   */
/* ---------------------------------------------------- */

int main(int argn, char *argv[]) {

    /* Declarations */

    int   i, pos;
    int   winmax, width, step, prior;
    int   carg, errflag, aflag, vflag; 

    MarkovMatrix *pmat;
    ProbaArray   win, baye;
    
    FastaSequence *seq;
    
    FILE  *filin;
    
    char  inname[FILENAME_MAX];

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
 
    aflag  = 0;         /* binary input         */
    vflag  = 0;         /* not verbose          */
    
    prior  = PCODANT;   /* 80% coding           */
    width  = WIN_WIDTH; /* window width         */
    step   = WIN_STEP;  /* window step          */

    (void) strcpy(inname, MATNAME);

    while ((carg = getopt(argn, argv, "ahm:p:s:vw:")) != -1) {
            
        switch(carg) {

            case 'a' :              /* ascii mode   */
                aflag = 1;
                break;

            case 'h' :              /* help         */
                sPrintHelp();
                exit(0);
                break;

            case 'm' :              /* matrix name  */
                (void) strcpy(inname, optarg);
                break;

            case 'p' :              /* prior        */
                if (sscanf(optarg, "%d", &prior) != 1)
                   errflag++;
                break;

            case 's' :              /* step         */
                if (sscanf(optarg, "%d", &step) != 1)
                   errflag++;
                break;

            case 'v' :              /* verbose      */
                vflag = 1;
                break;

            case 'w' :              /* width        */
                if (sscanf(optarg, "%d", &width) != 1)
                   errflag++;
                break;

            case '?' :              /* misusage     */
                errflag++;
        }
    }
    
    /* ------------------------ */
    /* may remain 1 argument    */

    if ((argn - optind) > 1)
        errflag++;
    else if ((argn - optind) == 1) {
        if (! AssignToStdin(argv[optind]))
           errflag++;
    }

    /* ------------------------ */
    /* check arguments          */
    
    if ((prior <= 0) || (prior >= 100)) {
        Erreur("invalid prior value", 0);
        errflag++;
    }

    if (width <= 0) {
        Erreur("invalid width value", 0);
        errflag++;
    }

    if (step <= 0) {
        Erreur("invalid step value", 0);
        errflag++;
    }

    if (! (filin = OpenFile(MatrixPathName(inname), (aflag ? "r" : "rb"))))
        errflag++;

    /* ------------------------ */
    /* exit on usage error      */

    if (errflag) 
        (void) sExitUsage(1);

    /* -------------------------------- */
    /* lecture matrice                  */
    /* -------------------------------- */

    if (! ReadMatrix(pmat, aflag, filin))
       exit(3);

    (void) fclose(filin);

    /* -------------------------------- */
    /* calcul des proba de transition   */
    /* et probas initiales              */
    /* -------------------------------- */
    
    ProbaCond(pmat);   

    /* -------------------------------- */
    /* lecture sequence                 */
    /* -------------------------------- */

    seq = NewFastaSequence();
    ReadFastaSequence(stdin, seq);
    (void) RNAtoDNA(Upper(seq->seq));

    /* -------------------------------- */
    /* boucle deplacement fenetre       */
    /* et ecriture des probas dans les  */
    /* 7 phases                         */
    /* -------------------------------- */

    winmax = width - MAX(pmat->kupleC, pmat->kupleN) + 1;

    for (pos = 0 ; pos < seq->length - width + 1 ; pos += step) {

       ProbaMarkov(seq->seq + pos, winmax, pmat, &win);

       ProbaBayes(&win, (float) prior / 100., &baye);

       (void) printf("%6d\t", pos + width/2 + 1);
       
       for (i = 0 ; i < 3 ; i++)
          (void) printf("%.3f\t", baye.probpos[InternalFrame(i, 1, seq->length)]);

       for (i = 0 ; i < 3 ; i++)
          (void) printf("%.3f\t", baye.probneg[InternalFrame(i, 0, seq->length)]);

       (void) printf("%.3f\n", baye.probnc);
    }

    /* -------------------------------- */
    /* fin du programme                 */
    /* -------------------------------- */

    FREE(pmat);

    FreeFastaSequence(seq);
    
    return 0 ;
}

