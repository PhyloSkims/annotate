/* ---------------------------------------------------------------- */
/* Copyright (c) Atelier de BioInformatique                         */
/* @file: prokov_score.c                                            */
/* @desc: prokov - phase detection : scoring                        */
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
    PP  "  * detection phase - score CDS's *            \n");
    PP  "                                               \n");
    PP  "usage: prokov_score [options] fastafile        \n");
    PP  "------------------------------------------     \n");
    PP  "options:                                       \n");
    PP  "                                               \n");      
    PP  "-a : read matrix in [A]scii format             \n");
    PP  "        default: read in binary format         \n");
    PP  "                                               \n");      
    PP  "-h : this [H]elp                               \n");
    PP  "                                               \n");      
    PP  "-I : [I]gnore last codon of each sequence      \n");
    PP  "        default: off                           \n");
    PP  "                                               \n");      
    PP  "-m name : [M]atrix filename                    \n");
    PP  "        default: -m %s\n", MATNAME);
    PP  "                                               \n");      
    PP  "-o : [O]utput sequence data too (fasta format) \n");
    PP  "        default: off (short output)            \n");
    PP  "                                               \n");      
    PP  "-p nn : use nn %% as prior coding probability  \n");
    PP  "        default: -p %d\n", PCODANT);
    PP  "                                               \n");      
    PP  "-S string : use string as [S]tops              \n");
    PP  "        string has the form \"/XYZ/.../XYZ/\"  \n");
    PP  "        default: -S %s\n", STOP_DFT);
    PP  "        note: this is only useful with         \n");
    PP  "              the -z option                    \n");
    PP  "                                               \n");      
    PP  "-t nn : lower [T]hreshold - output only        \n");
    PP  "        if P[frame1] >= nn %%                  \n");
    PP  "        default: -t 0                          \n");
    PP  "        note: may be combined with -T          \n");
    PP  "                                               \n");      
    PP  "-T nn : upper [T]hreshold - output only        \n");
    PP  "        if P[frame1] <= nn %%                  \n");
    PP  "        default: -T 100                        \n");
    PP  "        note: may be combined with -t          \n");
    PP  "                                               \n");      
    PP  "-v : set [V]erbose mode                        \n");
    PP  "                                               \n");      
    PP  "-z : force proba to -1 if stop codon           \n");
    PP  "        default: off                           \n");
    PP  "                                               \n");      
    PP  "------------------------------------------     \n");
    PP  " seqfile : file containing CDS sequences       \n");
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

    int   i, lenseq, prior, upthres, lothres;
    int   carg, errflag, aflag, vflag, zflag, iflag, oflag; 
    float zupthres, zlothres;

    MarkovMatrix *pmat;
    ProbaArray   win, baye;
    
    FastaSequence *seq;
    
    FILE  *filin;
    
    char  stops[256], cstops[256], buffer[BUFSIZ], bufnum[BUFSIZ];

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
    iflag  = 0;         /* use last codon       */
    vflag  = 0;         /* not verbose          */
    zflag  = 0;         /* no zero stops        */
    oflag  = 0;         /* short output         */
    
    prior  = PCODANT;   /* 80% coding           */
    upthres = 100;      /* upper threshold      */
    lothres = 0;        /* lower threshold      */

    (void) strcpy(inname, MATNAME);
    (void) strcpy(stops, STOP_DFT);

    while ((carg = getopt(argn, argv, "ahIm:op:S:t:T:vz")) != -1) {
            
        switch(carg) {

            case 'a' :              /* ascii mode   */
                aflag = 1;
                break;

            case 'h' :              /* help         */
                sPrintHelp();
                exit(0);
                break;

            case 'I' :              /* Ignore last  */
                iflag = 1;
                break;

            case 'm' :              /* matrix name  */
                (void) strcpy(inname, optarg);
                break;

            case 'o' :              /* long output */
                oflag = 1;
                break;

            case 'p' :              /* prior        */
                if (sscanf(optarg, "%d", &prior) != 1)
                   errflag++;
                break;

            case 'S' :              /* stops        */
                (void) strcpy(stops, optarg);
                break;

            case 't' :              /* low thresh.  */
                if (sscanf(optarg, "%d", & lothres) != 1)
                   errflag++;
                break;

            case 'T' :              /* up thresh.   */
                if (sscanf(optarg, "%d", & upthres) != 1)
                   errflag++;
                break;

            case 'v' :              /* verbose      */
                vflag = 1;
                break;

            case 'z' :              /* zero stop    */
                zflag = 1;
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

    if ((lothres < 0) || (lothres > 100)) {
        Erreur("invalid lower threshold (-o) value", 0);
        errflag++;
    }

    if ((upthres < 0) || (upthres > 100)) {
        Erreur("invalid upper threshold (-O) value", 0);
        errflag++;
    }

    if (upthres < lothres) {
        Erreur("incoherent lower/upper threshold (-o -O) values", 0);
        errflag++;
    }

    if (! (filin = OpenFile(MatrixPathName(inname), (aflag ? "r" : "rb"))))
        errflag++;

    /* ------------------------ */
    /* exit on usage error      */

    if (errflag) 
        (void) sExitUsage(1);

    /* ------------------------ */
    /* initialize codons        */
    
    Complement(Reverse(strcpy(cstops,  stops)));

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

    zlothres = ((float) lothres) / 100.;
    zupthres = ((float) upthres) / 100.;

    /* -------------------------------- */
    /* boucle lecture des sequences     */
    /* -------------------------------- */

    seq = NewFastaSequence();

    while (ReadFastaSequence(stdin, seq)) {

       if (! seq->ok) {
           (void) sprintf(buffer, "invalid sequence %s", seq->name);
           Erreur(buffer, 0);
           continue;
       }
       
       if (vflag)
          (void) fprintf(stderr, "%s [%ld b]\n", seq->name, seq->length);
        
       RNAtoDNA(Upper(seq->seq));

       if (iflag && (seq->length >= 4)) {
          seq->seq[seq->length-4] = '\000';
          seq->length -= 3;
       }
       
       lenseq = seq->length - MAX(pmat->kupleC, pmat->kupleN) + 1;
       
       ProbaMarkov(seq->seq, lenseq, pmat, &win);

       ProbaBayes(&win, (float) prior / 100., &baye);

       if (zflag) {                             /* -1 score when stop */

          for (i = 0 ; i < 3 ; i++) {
             if (FindCodon(seq->seq, seq->length, stops, i, seq->length - 1) >= 0)
                baye.probpos[i] = -1.;
          }       
          
          for (i = 0 ; i < 3 ; i++) {
             if (FindCodon(seq->seq, seq->length, cstops, i, seq->length - 1) >= 0)
                baye.probneg[i] = -1.;
          }
       }
       
                                                /* output sequence      */

       if ((baye.probpos[0] >= zlothres) && (baye.probpos[0] <= zupthres)) {

           (void) strcpy(buffer, " | ");

           for (i = 0 ; i < 3 ; i++) {
              (void) sprintf(bufnum, "%6.3f ", baye.probpos[InternalFrame(i, 1, seq->length)]);
              (void) strcat(buffer, bufnum);
           }
           for (i = 0 ; i < 3 ; i++) {
              (void) sprintf(bufnum, "%6.3f ", baye.probneg[InternalFrame(i, 0, seq->length)]);
              (void) strcat(buffer, bufnum);
           }
           (void) sprintf(bufnum, "%6.3f", baye.probnc);
           (void) strcat(buffer, bufnum);

           if (oflag) {
              seq->comment[MAX(0, FASTA_COMLEN - strlen(buffer) - 1)] = '\000';
              (void) strcat(seq->comment, buffer);         
              WriteFastaSequence(stdout, seq, FASTA_CHAR_PER_LINE);
           }
           else
                printf("%s %s\n", seq->name, buffer + 3);
       }

    }

    /* -------------------------------- */
    /* fin du programme                 */
    /* -------------------------------- */

    FREE(pmat);

    FreeFastaSequence(seq);
    
    return 0 ;
}

