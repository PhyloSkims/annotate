/* ---------------------------------------------------------------- */
/* Copyright (c) Atelier de BioInformatique                         */
/* @file: prokov_orf.c                                              */
/* @desc: prokov - : recherche des orfs                             */
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

#define BIOFRAME 1

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
    PP  "  * simple orfs searching utility *            \n");
    PP  "                                               \n");
    PP  "usage: prokov_orf [options] seqfile            \n");
    PP  "------------------------------------------     \n");
    PP  "options:                                       \n");
    PP  "                                               \n");      
    PP  "-h : this [H]elp                               \n");
    PP  "                                               \n");      
    PP  "-l nn : minimun cds [L]ength                   \n");
    PP  "        default: -l %d\n", MIN_CDS);
    PP  "                                               \n");      
    PP  "-o : [O]utput sequence data too (fasta format) \n");
    PP  "        default: off (short output)            \n");
    PP  "                                               \n");      
    PP  "-s string : use string as Starts               \n");
    PP  "        string has the form \"/XYZ/.../XYZ/\"  \n");
    PP  "        default: -S %s\n", START_DFT);
    PP  "                                               \n");      
    PP  "-S string : use string as Stops                \n");
    PP  "        string has the form \"/XYZ/.../XYZ/\"  \n");
    PP  "        default: -S %s\n", STOP_DFT);
    PP  "                                               \n");      
    PP  "-v : set [V]erbose mode                        \n");
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
    PP  "usage: prokov_orf [options] seqfile    \n");
    PP  "type \"prokov_orf -h\" for help        \n");

    if (stat)
        exit(stat);
            
    return stat;
}

#undef PP

/* ----------------------------------------------- */
/* output fragment of sequence (fasta format)      */
/* ----------------------------------------------- */
void sOutputSequence(FastaSequence *seq, int posa,
                     int poss, int strand, 
                     char *name, char *comment)
{
    char save;
        
    static FastaSequence sSeq;

    (void) strcpy(sSeq.name, name);
    (void) strcpy(sSeq.comment, comment);

    sSeq.seq = seq->seq + posa;

    poss += 3;

    save = seq->seq[poss];
    seq->seq[poss] = '\000';
    sSeq.length = poss - posa;

    WriteFastaSequence(stdout, &sSeq, FASTA_CHAR_PER_LINE);

    seq->seq[poss] = save;
}

/* ----------------------------------------------- */
/* print a CDS                                     */
/* ----------------------------------------------- */

void sPrintCds(FastaSequence *seq, int strand, int phase,
               int posa, int poss, int outflag)
{
    int from, to, frame;
    char bufname[BUFSIZ], bufcom[BUFSIZ];
    
    from = ((strand > 0) ? posa + 1 : seq->length - poss - 2);
    to   = ((strand > 0) ? poss + 3 : seq->length - posa);

#if BIOFRAME
    frame = ((from - 1) % 3 ) + 1;
#else
    frame = phase;
#endif

    (void) sprintf(bufname, "%s_ORF_%d_%d_%s", seq->name,
                               from, to, ((strand > 0) ? "D" : "C"));  

    (void) sprintf(bufcom, "%d %d %1s%d",
                   from, to, ((strand > 0) ? "+" : "-"), frame);

     if (outflag)       /* long output */    
        sOutputSequence(seq, posa, poss, strand, bufname, bufcom);
     else               /* short output */
        (void) printf("%s %s\n", bufname, bufcom);
}


/* ----------------------------------------------- */
/* process cds's in single phase                   */
/* ----------------------------------------------- */
static void sProcessSequence(FastaSequence *seq, 
                int strand, int phase, int mincds, int outflag,
                char *starts, char *stops)
{
    int posa, poss, lencds;


    posa = phase;                       /* start on phase */
    
    while ((posa = FindCodon(seq->seq, seq->length, starts, posa, seq->length - 1)) >= 0)  {

        poss = FindCodon(seq->seq, seq->length, stops, posa, seq->length - 1);

        if (poss < 0)
           break;               /* no more stop */
                
        lencds = poss - posa;
             
        if (lencds < mincds) {
           posa = poss + 3;     /* try next start */
           continue;            
        }
        
        sPrintCds(seq, strand, phase, posa, poss, outflag); 
        posa = poss + 3;

    }
}

/* ---------------------------------------------------- */
/*                PROGRAMME PRINCIPAL                   */
/* ---------------------------------------------------- */

int main(int argn, char *argv[]) {

    /* Declarations */

    int   phase, mincds;
    int   carg, errflag, vflag, oflag; 

    FastaSequence *seq;
    
    FILE  *filin;

    char  stops[256], starts[256];
    
    /* ------------------------ */
    /* lecture des arguments    */
    /* ------------------------ */

    /* ------------------------ */
    /* valeurs par defaut       */
    
    errflag = 0;
 
    vflag  = 0;         /* not verbose          */
    oflag  = 0;         /* short output         */
    
    mincds = MIN_CDS;   /* min cds length       */
    
    (void) strcpy(starts, START_DFT);
    (void) strcpy(stops,  STOP_DFT);

    while ((carg = getopt(argn, argv, "hl:os:S:v")) != -1) {
            
        switch(carg) {

            case 'h' :              /* help         */
                sPrintHelp();
                exit(0);
                break;

            case 'l' :              /* min cds      */
                if (sscanf(optarg, "%d", &mincds) != 1)
                   errflag++;
                break;

            case 'o' :              /* verbose      */
                oflag = 1;
                break;

            case 's' :              /* starts       */
                (void) strcpy(starts, optarg);
                break;

            case 'S' :              /* stops        */
                (void) strcpy(stops, optarg);
                break;

            case 'v' :              /* verbose      */
                vflag = 1;
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
    
    if (mincds <= 0) {
        Erreur("invalid mincds value", 0);
        errflag++;
    }

    /* ------------------------ */
    /* exit on usage error      */

    if (errflag) 
        (void) sExitUsage(1);

    /* -------------------------------- */
    /* lecture sequence                 */
    /* -------------------------------- */

    seq = NewFastaSequence();
    ReadFastaSequence(stdin, seq);
    (void) RNAtoDNA(Upper(seq->seq));

    /* -------------------------------- */
    /* recherche des cds dans chaque    */
    /* phase positive                   */
    /* -------------------------------- */
    
    for (phase = 0 ; phase < 3 ; phase++)
        sProcessSequence(seq, 1, phase, mincds, oflag, starts, stops);

    /* -------------------------------- */
    /* recherche des cds dans chaque    */
    /* phase negative                   */
    /* -------------------------------- */

    Complement(Reverse(seq->seq));

    for (phase = 0 ; phase < 3 ; phase++)
        sProcessSequence(seq, -1, phase, mincds, oflag, starts, stops);

    /* -------------------------------- */
    /* fin du programme                 */
    /* -------------------------------- */

    FreeFastaSequence(seq);
    
    return 0 ;
}

