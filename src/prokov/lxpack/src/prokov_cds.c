/* ---------------------------------------------------------------- */
/* Copyright (c) Atelier de BioInformatique                         */
/* @file: prokov_cds.c                                              */
/* @desc: prokov - phase detection : recherche cds                  */
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
    PP  "  * detection phase - CDSs searching *         \n");
    PP  "                                               \n");
    PP  "usage: prokov_cds [options] seqfile            \n");
    PP  "------------------------------------------     \n");
    PP  "options:                                       \n");
    PP  "                                               \n");      
    PP  "-a : read matrix in [A]scii format             \n");
    PP  "        default: read in binary format         \n");
    PP  "                                               \n");      
    PP  "-h : this [H]elp                               \n");
    PP  "                                               \n");      
    PP  "-l nn : minimun cds [L]ength                   \n");
    PP  "        default: -l %d\n", MIN_CDS);
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
    PP  "-s string : use string as Starts               \n");
    PP  "        string has the form \"/XYZ/.../XYZ/\"  \n");
    PP  "        default: -S %s\n", START_DFT);
    PP  "                                               \n");      
    PP  "-S string : use string as Stops                \n");
    PP  "        string has the form \"/XYZ/.../XYZ/\"  \n");
    PP  "        default: -S %s\n", STOP_DFT);
    PP  "                                               \n");      
    PP  "-t nn : cds probability [T]hreshold  (%%)      \n");
    PP  "        default: -l %d\n", THRESH_CDS);
    PP  "                                               \n");      
    PP  "-v : set [V]erbose mode                        \n");
    PP  "                                               \n");      
    PP  "-w nn : start [W]alk                           \n");
    PP  "        expressed as %% of cds size            \n");           
    PP  "        default: walk off                      \n");
    PP  "                                               \n");      
    PP  "-W nn : start [W]alk                           \n");
    PP  "        expressed as number of bases           \n");           
    PP  "        default: walk off                      \n");
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
    PP  "usage: prokov_cds [options] seqfile    \n");
    PP  "type \"prokov_cds -h\" for help        \n");

    if (stat)
        exit(stat);
            
    return stat;
}

#undef PP

/* ----------------------------------------------- */
/* score a single window in phase +1               */
/* ----------------------------------------------- */
static float sScoreWindow(char *seq, int lenseq, MarkovMatrix *mat,
                          int prior)
{
    ProbaArray   win, bay;

    ProbaMarkov(seq, lenseq, mat, &win);
    ProbaBayes(&win, (float) prior / 100., &bay);
    
    return bay.probpos[0];      /* phase +1 */
}

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
               int posa, int poss, int poso, float score,
               int outflag)
{
    int from, to, ori, frame;
    char bufname[BUFSIZ], bufcom[BUFSIZ];
    

    from = ((strand > 0) ? posa + 1 : seq->length - poss - 2);
    to   = ((strand > 0) ? poss + 3 : seq->length - posa);

    ori  = ((strand > 0) ? poso + 1 : seq->length - poso);

#if BIOFRAME
    frame = ((from - 1) % 3 ) + 1;
#else
    frame = phase;
#endif

    (void) sprintf(bufname, "%s_CDS_%d_%d_%s", seq->name,
                               from, to, ((strand > 0) ? "D" : "C"));  

    (void) sprintf(bufcom, "%d %d %1s%d %.3f %s %d",
                   from, to, 
                   ((strand > 0) ? "+" : "-"), frame,
                   score,
                   ((poso == posa) ? "max" : "walk"),
                   ((poso == posa) ? 0 : ori));

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
                int dowalk, int walk, int prior, int thresh,
                char *starts, char *stops,
                MarkovMatrix *mat)
{
    int posa, posb, poso, poss, lencds, awalk;
    float scort, scora, scorb;


    scort = (float) thresh / 100.;      /* threshold      */
    
    posa = phase;                       /* start on phase */
    
    awalk = (  dowalk 
             ? ((walk > 0) ? walk : ABS(walk) * seq->length / 100)
             : 0);

    while ((posa = FindCodon(seq->seq, seq->length, starts, posa, seq->length - 1)) >= 0)  {

        poss = FindCodon(seq->seq, seq->length, stops, posa, seq->length - 1);

        if (poss < 0)
           break;               /* no more stop */
                
        lencds = poss - posa;
             
        if (lencds < mincds) {
           posa = poss + 3;     /* try next start */
           continue;            
        }
        
        scora = sScoreWindow(seq->seq + posa, lencds, mat, prior);

        /* try to walk around start */
        
        poso = posa;
        
        if (dowalk) {
        
            posb = posa + 3;
            
            while ((posb = FindCodon(seq->seq, seq->length, starts, posb, seq->length - 1)) >= 0) {
            
                if ((posb - posa) > awalk)
                   break;       /* too far away */

                lencds = poss - posb;
            
                if (lencds < mincds)
                   break;       /* too short    */
            
                scorb = sScoreWindow(seq->seq + posb, lencds, mat, prior);
        
                if (scorb > scora) {  /* best start */
                   posa = posb;
                   scora = scorb;
                }
                
                posb += 3;

             }
        }

        if (scora >= scort) {
           sPrintCds(seq, strand, phase, posa, poss, poso, scora, outflag); 
           posa = poss + 3;
        }
        else {
           posa += 3;
        }
    }
}

/* ---------------------------------------------------- */
/*                PROGRAMME PRINCIPAL                   */
/* ---------------------------------------------------- */

int main(int argn, char *argv[]) {

    /* Declarations */

    int   phase, walk, mincds, prior, thresh;
    int   carg, errflag, aflag, vflag, wflag, oflag; 

    MarkovMatrix *pmat;
  
    FastaSequence *seq;
    
    FILE  *filin;

    char  stops[256], starts[256];
    
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
    wflag  = 0;         /* no walk              */
    oflag  = 0;         /* short output         */
    
    prior  = PCODANT;   /* 80% coding           */
    walk   = 0;         /* walk off             */
    mincds = MIN_CDS;   /* min cds length       */
    thresh = THRESH_CDS;/* 40                   */
    
    (void) strcpy(inname, MATNAME);
    (void) strcpy(starts, START_DFT);
    (void) strcpy(stops,  STOP_DFT);

    while ((carg = getopt(argn, argv, "ahl:m:op:s:S:t:vw:W:")) != -1) {
            
        switch(carg) {

            case 'a' :              /* ascii mode   */
                aflag = 1;
                break;

            case 'h' :              /* help         */
                sPrintHelp();
                exit(0);
                break;

            case 'l' :              /* min cds      */
                if (sscanf(optarg, "%d", &mincds) != 1)
                   errflag++;
                break;

            case 'm' :              /* matrix name  */
                (void) strcpy(inname, optarg);
                break;

            case 'o' :              /* long output    */
                oflag = 1;
                break;

            case 'p' :              /* prior        */
                if (sscanf(optarg, "%d", &prior) != 1)
                   errflag++;
                break;

            case 's' :              /* starts       */
                (void) strcpy(starts, optarg);
                break;

            case 'S' :              /* stops        */
                (void) strcpy(stops, optarg);
                break;

            case 't' :              /* threshold    */
                if (sscanf(optarg, "%d", &thresh) != 1)
                   errflag++;
                break;

            case 'v' :              /* verbose      */
                vflag = 1;
                break;

            case 'w' :              /* walk         */
                if (sscanf(optarg, "%d", &walk) != 1)
                   errflag++;
                wflag = 1;
                walk = -walk;      /* means '%'     */
                break;

            case 'W' :              /* walk         */
                if (sscanf(optarg, "%d", &walk) != 1)
                   errflag++;
                wflag = 1;
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

    if ((thresh < 0) || (thresh > 100)) {
        Erreur("invalid thresh value", 0);
        errflag++;
    }

    if (mincds <= 0) {
        Erreur("invalid mincds value", 0);
        errflag++;
    }

    if ((wflag) && (walk < 0) && (ABS(walk) > 100)) {
       Erreur("invalid walk value", 0);
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
    /* recherche des cds dans chaque    */
    /* phase positive                   */
    /* -------------------------------- */
    
    for (phase = 0 ; phase < 3 ; phase++)
        sProcessSequence(seq, 1, phase, mincds, oflag, wflag,
                         walk, prior, thresh, starts, stops,
                         pmat);

    /* -------------------------------- */
    /* recherche des cds dans chaque    */
    /* phase negative                   */
    /* -------------------------------- */

    Complement(Reverse(seq->seq));

    for (phase = 0 ; phase < 3 ; phase++)
        sProcessSequence(seq, -1, phase, mincds, oflag, wflag,
                         walk, prior, thresh, starts, stops,
                         pmat);

    /* -------------------------------- */
    /* fin du programme                 */
    /* -------------------------------- */

    FREE(pmat);

    FreeFastaSequence(seq);
    
    return 0 ;
}

