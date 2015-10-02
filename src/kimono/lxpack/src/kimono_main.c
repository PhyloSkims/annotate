/* ==================================================== */
/* @File kimono_main.c                                  */
/* @history:                                            */
/* @+   Apr. 97 <Gloup> first draft                     */
/* ==================================================== */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "Genetic.h"  
#include "fasta_io.h"
#include "kimono.h"

#ifdef HAS_UNISTD_H
#include HAS_UNISTD_H
#endif

#ifdef HAS_GETOPT_H
#include HAS_GETOPT_H
#endif

#ifndef FILENAME_MAX
#define FILENAME_MAX    1024    /* max # of characters in a path name */
#endif

#ifndef BUFSIZ
#define BUFSIZ          1024
#endif

#define SELECTOR_NONE   0
#define SELECTOR_COUNT  1
#define SELECTOR_CUMUL  2
#define SELECTOR_CODON  3
#define SELECTOR_GENE   4

#define DISTANCE_NONE   0
#define DISTANCE_CHI    1
#define DISTANCE_EUCLID 2
#define DISTANCE_HELL   3

#define MIXER_NONE      0
#define MIXER_BASIC     1
#define MIXER_CROSS     2

#define DFT_WINDOW      -25  /* aka 25% */
#define DFT_STEP        1000

#define CUMUL_NORM      1.0e6

#define DEBUG           0

/* -------------------------------------------- */
/* macros                                       */
/* -------------------------------------------- */

/* -------------------------------------------- */
/* globals                                      */
/* -------------------------------------------- */

static char sDNA[]  = KIM_DNA;

static float sCodonScore[NB_NUC][NB_CODONS];

/* -------------------------------------------- */
/* getopt globals                               */
/* -------------------------------------------- */

extern char *optarg;
extern int optind;

/* =============================================== */
/* Utilities Functions   			   */
/* =============================================== */

/* ----------------------------------------------- */
/* make a count selector                           */
/* ----------------------------------------------- */

static int sSetSelector(char *buf) {

  if (! strcmp(buf, "count"))
    return SELECTOR_COUNT;

  if (! strcmp(buf, "cumul"))
    return SELECTOR_CUMUL;

  if (! strcmp(buf, "codon"))
    return SELECTOR_CODON;

  if (! strcmp(buf, "genes"))
    return SELECTOR_GENE;
    
  return SELECTOR_NONE;
}

/* ----------------------------------------------- */
/* get the count selector                          */
/* ----------------------------------------------- */

static char *sGetSelector(int selector) {

  switch (selector) {
  
    case SELECTOR_COUNT : return "count";
    case SELECTOR_CUMUL : return "cumul";
    case SELECTOR_CODON : return "codon";
    case SELECTOR_GENE  : return "genes";

  }
  
  return "<unknown>";
}

/* ----------------------------------------------- */
/* make a distance selector                        */
/* ----------------------------------------------- */

static int sSetDistance(char *buf) {

  if (! strcmp(buf, "chi"))
    return DISTANCE_CHI;

  if (! strcmp(buf, "euclid"))
    return DISTANCE_EUCLID;

  if (! strcmp(buf, "hell"))
    return DISTANCE_HELL;

  return DISTANCE_NONE;
}

/* ----------------------------------------------- */
/* get the count selector                          */
/* ----------------------------------------------- */

static char *sGetDistance(int selector) {

  switch (selector) {
  
    case DISTANCE_CHI    : return "chi";
    case DISTANCE_EUCLID : return "euclid";
    case DISTANCE_HELL   : return "hell";
  }
  
  return "<unknown>";
}

/* ----------------------------------------------- */
/* make a mixer selector                           */
/* ----------------------------------------------- */

static int sSetMixer(char *buf) {

  if (! strcmp(buf, "basic"))
    return MIXER_BASIC;

  if (! strcmp(buf, "cross"))
    return MIXER_CROSS;

  return MIXER_NONE;
}

/* ----------------------------------------------- */
/* get the mixer selector                          */
/* ----------------------------------------------- */

static char *sGetMixer(int selector) {

  switch (selector) {
  
    case MIXER_BASIC     : return "basic";
    case MIXER_CROSS     : return "cross";
  }
  
  return "<unknown>";
}

/* ----------------------------------------------- */
/* make a phase mask                               */
/* ----------------------------------------------- */

static int sGetPhaseMask(char *buf) {
 
  int mask = 0x0;
  
  if (strchr(buf, '1'))
    mask |= 0x1;

  if (strchr(buf, '2'))
    mask |= 0x2;

  if (strchr(buf, '3'))
    mask |= 0x4;

  return mask;
}

/* ----------------------------------------------- */
/* uppercase string                                */
/* ----------------------------------------------- */

static char *sUpper(char *buf) {

  char *c;

  for (c = buf ; *c ; c++) {
    *c = UPPER(*c);
  }
  
  return buf;
}

/* ----------------------------------------------- */
/* replace symbols in alphaOut by MASKED_SYMBOL    */
/* ----------------------------------------------- */

static char *sExcludeFromSeq(char *buf, char *alphaOut, char maskSymbol) {

  char *c;
  
  if (! *alphaOut)
    return buf;

  for (c = buf ; *c ; c++) {
    if (strchr(alphaOut, *c))
      *c = maskSymbol;
  }
  
  return buf;
}

/* ----------------------------------------------- */
/* check alphabet                                  */
/* ----------------------------------------------- */

static int sCheckAlphabet(char *alpha) {

  char *c;
  for (c = alpha ; *c ; c++) {
    if (! strchr(sDNA, *c))   /* invalid symbol   */
      return 0;
    if (strchr(c+1, *c))      /* duplicate symbol */
      return 0;
  }
  
  return strlen(alpha);
}

/* ----------------------------------------------- */
/* main                                            */
/* ----------------------------------------------- */

main(int argn, char *argv[])
{
  int            i, j, carg, errflag;
  int            gflag, iflag, xflag, dflag, code;
  int            selector, mixer, phase;
  int            pos, len, step, window, nopen, nvalid, nbins, sign;
  int            allLeft[3 * NB_NUC + 1][2], allRight[3 * NB_NUC + 1][2],
                 countLeft[3 * NB_NUC + 1], countRight[3 * NB_NUC + 1],
                 tmpLeft[2], tmpRight[2];
  float          skew, sum;

  Alphabet       *alpha;
  Sequence       *seq;
  SetOfIntervals *zones;
  char           alphaOut[128];
  char           inputfile[FILENAME_MAX], includefile[FILENAME_MAX];

  
  IndexFunction  indexingFunction;
  CountFunction  countingFunction;
  DistFunction   distanceFunction;
  MixFunction    mixingFunction;
  

  /* ---------------------------- */    
  /* defaut options               */

  errflag     = 0;
  selector    = SELECTOR_COUNT;  /* mono countings      */
  mixer       = MIXER_BASIC;     /* basic mixer         */
  dflag       = DISTANCE_CHI;    /* chi-square distance */
  iflag       = 0;               /* no inclusion        */
  phase       = 0x7;             /* all phases          */
  xflag       = 0;               /* no exclusion        */
  code        = 0;               /* universal code      */
  step        = DFT_STEP;
  window      = DFT_WINDOW;

  alpha = NucleicAlphabet();
  
  alphaOut[0] = '\000';
  
  (void) strcpy(inputfile, "<stdin>");
  (void) strcpy(includefile, "<none>");
  
  /* ---------------------------- */    
  /* parse cmdline arguments      */
        
  while ((carg = getopt(argn, argv, "a:b:c:d:e:g:hi:m:p:s:w:W:x:")) != -1) {
            
    switch(carg) {

      case 'a' :                        /* [a]lphabet     */
        (void) FreeAlphabet(alpha);
        alpha = NewAlphabet(sUpper(optarg));
        if (! sCheckAlphabet(alpha->symbol)) {
          Erreur("bad alphabet", -1);
          errflag++;
        }
        break;

      case 'b' :                        /* excluded symbols */
        (void) strcpy(alphaOut, sUpper(optarg));
        break;

    
      case 'c' :                        /* [c]ount selector */
        selector = sSetSelector(optarg);
        if (selector == SELECTOR_NONE) {
          Erreur("bad count selector", -1);
          errflag++;
        }
        break;

      case 'e':                           /* g[e]netic code   */
        if (   (sscanf(optarg, "%d", &code) != 1)
            || (code < 0) || (code > GEN_MAX_CODES)) {
          Erreur("bad code value: -e (0-8)", -1);
          errflag++;
        }
        else {
          InitCodonSkewCounter(code);
        }
        break;

      case 'd' :                        /* [d]istance selector */
        dflag = sSetDistance(optarg);
        if (dflag == DISTANCE_NONE) {
          Erreur("bad distance selector", -1);
          errflag++;
        }
        break;

      case 'h' :                        /* help             */
        PrintHelp();
        exit(0);

      case 'i' :                        /* [i]nclude file   */
        if (xflag) {
          Erreur("-x and -i options are incompatibles", -1);
          errflag++;
        }
        iflag = 1;
        (void) strcpy(includefile, optarg);
        break;

      case 'm' :                        /* [m]ixer selector */
        mixer = sSetMixer(optarg);
        if (mixer == MIXER_NONE) {
          Erreur("bad mixer selector", -1);
          errflag++;
        }
        break;

      case 'p' :                        /* [p]hase mask     */
        phase = sGetPhaseMask(optarg);
        if (phase == 0) {
          Erreur("bad phase value", -1);
          errflag++;
        }
        break;

      case 's':                         /* [s]tep           */
        if (   (sscanf(optarg, "%d", &step) != 1)
            || (step <= 0) ) {
          Erreur("bad step value", 0);
          errflag++;
        }
        break;

      case 'w':                         /* [w]indow         */
        if (   (sscanf(optarg, "%d", &window) != 1)
            || (window <= 0) ) {
          Erreur("bad window value", -1);
          errflag++;
        }
        break;

      case 'W':                        /* [W]indow          */
        if ( (sscanf(optarg, "%d", &window) != 1)
            || (window <= 0) 
            || (window > 100)) {
          Erreur("bad bad window percent value", -1);
          errflag++;
        }
        window = -window;
        break;

      case 'x' :                        /* e[X]lude file   */
        if (iflag) {
          Erreur("-x and -i options are incompatibles", -1);
          errflag++;
        }
        xflag = 1;
        (void) strcpy(includefile, optarg);
        break;
              
      case '?' :                        /* misusage         */
        errflag++;

    }
  }

  /* -------------------------------- */    
  /* check some options compatibility */
  
  if ((phase != 0x7) && (! iflag)) {
    Erreur("-p flag requires -i flag", -1);
    errflag++;
  }

  if (selector == SELECTOR_CUMUL) {
    if (alpha->length != 2) {
      Erreur("need a two-letters alphabet for -c cumul", -1);
      errflag++;
    }
  }

  if (selector == SELECTOR_GENE) {
    if (! iflag) {
      Erreur("-c gene requires -i flag", -1); 
      errflag++;
    }
    if (phase != 0x7) {
      Erreur("-c gene needs -p 123", -1);
      errflag++;
    }
    if (strcmp(alpha->symbol, sDNA)) {
      Erreur("-a ignored for -c gene", 0);
      (void) FreeAlphabet(alpha);
      alpha = NucleicAlphabet();
    }
  }

  if (selector == SELECTOR_CODON) {
    if (! iflag) {
      Erreur("-c codon requires -i flag", -1); 
      errflag++;
    }
    if (phase != 0x7) {
      Erreur("-c codon is incompatible with -p flag", -1);
      errflag++;
    }
    if (strcmp(alpha->symbol, sDNA)) {
      Erreur("-a ignored for -c codon", 0);
      (void) FreeAlphabet(alpha);
      alpha = NucleicAlphabet();
    }
  }
  
  /* ---------------------------- */    
  /* may remain 1 argument        */

  argn -= optind;
  
  if (argn > 1)
    errflag++;

  if (errflag) 
      ExitUsage(1);

  /* ---------------------------- */    
  /* allocate sequences           */
          
  seq = NewSequence(alpha);

  /* ---------------------------- */    
  /* read genome                  */

  if (argn == 1) {
    (void) strcpy(inputfile, argv[optind]);
    if (! AssignToStdin(inputfile))
      Erreur("cannot open sequence file", 3);
  }
   
  if (! ReadSequence(stdin, seq))
    Erreur("cannot read sequence file", 4);
    
  (void) sUpper(seq->fasta->seq);    /* upper case sequence  */
  
  (void) sExcludeFromSeq(seq->fasta->seq, alphaOut, MASKED_SYMBOL);
    
  /* ---------------------------- */    
  /* read include file            */
    
  if (iflag || xflag) {
    if (! (AssignToStdin(includefile)))
      Erreur("cannot open include file", 4);
    if (! (zones = ReadZones(seq->fasta->length)))
      Erreur("cannot read include zones", 5);
  }
  else {
    zones = NULL;
  }
  
  /* ---------------------------- */    
  /* adjust parameters            */
  
                                            /* -------------------- */
  len = seq->fasta->length;                 /* sequence length      */
  
                                            /* -------------------- */
  if (window < 0)                           /* window               */
    window = len * (-window) / 100;

  window /= 2;						 /* actual wondow for each half */
  
                                            /* -------------------- */
  if (selector == SELECTOR_GENE) {          /* alphabet             */
    (void) strcpy(alpha->symbol, "DR");
    alpha->length = 2;
  }
                                            /* -------------------- */
  nbins = alpha->length;                    /* # bins (for chi2)    */
                                            /* -------------------- */
  if (xflag) {                              /* indexing function    */
    indexingFunction = ExcludeIndexer;
  }
  else if (iflag) {
    indexingFunction = PhaseIndexer;
  }
  else {
    indexingFunction = NoZoneIndexer;
  }
  
  if (selector == SELECTOR_GENE) {
    indexingFunction = GeneIndexer;  
  }
  if (selector == SELECTOR_CODON) {
    indexingFunction = CodonIndexer;
  }

                                            /* -------------------- */
  if (selector == SELECTOR_CODON) {         /* counting function    */

    /* precompute score for each possible codon */
    /* to speedup computation                   */
    
    CodonSkew(code, sCodonScore);

    countingFunction = CodonSkewCounter;
    
    nbins = 3 * NB_NUC;
  }
  else {
    countingFunction = BasicCounter;
  }

                                            /* -------------------- */
  if (dflag == DISTANCE_EUCLID) {           /* distance function    */
    distanceFunction = Euclid;
  }
  else if (dflag == DISTANCE_HELL) {
    distanceFunction = HellTwo;
  }
  else {
    distanceFunction = ChisTwo;
  }
  
  
                                            /* -------------------- */
  if (mixer == MIXER_CROSS) {               /* mixing function      */
    mixingFunction = CrossMixer;
  }
  else {
    mixingFunction = BasicMixer;
  }
  
  /* ---------------------------- */    
  /* index sequence               */

   indexingFunction(seq, zones, phase);

  /* ---------------------------- */    
  /* # included positions         */

  for (i = nopen = nvalid = 0 ; i < len ; i++) {     
    if ((seq->cell[i].mask & EXCLUDE_MASK) == 0) {
      nopen++;
      if (seq->cell[i].index >= 0)
        nvalid++;
    }
  }

    
  /* ---------------------------- */    
  /* printout header              */
  
  printf("# input     : %s\n", inputfile);
  printf("# %s   : %s\n", (xflag ? "exclude" : "include"), includefile);
  printf("# selector  : %s\n", sGetSelector(selector));
  printf("# mixer     : %s\n", sGetMixer(mixer));
  printf("# distance  : %s\n", sGetDistance(dflag));
  printf("# phase     : %d\n", phase);
  printf("# alphabet  : %s\n", alpha->symbol);
  printf("# excluded  : X%s\n", alphaOut);
  printf("# seqname   : %s\n", seq->fasta->name);
  printf("# seqlength : %d\n", len);
  printf("# unmasked  : %d\n", nopen);
  printf("# valid     : %d\n", nvalid);
  printf("# window    : %d\n", 2 * window);
  printf("# step      : %d\n", step);
  printf("#\n");
  
  printf("# pos");

  if (selector == SELECTOR_COUNT) {
    for (i = 0 ; i < alpha->length ; i++) {
      printf(" Left%c Right%c", alpha->symbol[i], alpha->symbol[i]);
    }
    printf(" LeftOther RightOther");
    printf(" %s%s", sGetDistance(dflag), alpha->symbol);
    printf(" %s[%s]Other", sGetDistance(dflag), alpha->symbol);
    for (i = 0 ; i < alpha->length ; i++) {
      printf(" %s%c", sGetDistance(dflag), alpha->symbol[i]);
    }
  }
  
  if (selector == SELECTOR_CUMUL) {
    printf(" %c-%c/%c+%c", alpha->symbol[0], alpha->symbol[1], alpha->symbol[0], alpha->symbol[1]);
  }

  if (selector == SELECTOR_GENE) {
    printf(" LeftDirect RightDirect LeftReverse RightReverse");
    printf(" %s_Direct", sGetDistance(dflag));
  }

  if (selector == SELECTOR_CODON) {
    printf(" LeftTotalCodon RightTotalCodon");
    for (i = 0 ; i < alpha->length ; i++) {
        printf(" Left%c0 Left%cP Left%cM", alpha->symbol[i], alpha->symbol[i], alpha->symbol[i]);
        printf(" Right%c0 Right%cP Right%cM", alpha->symbol[i], alpha->symbol[i], alpha->symbol[i]);
    }
    for (i = 0 ; i < alpha->length ; i++) {
      printf(" Total%s%c", sGetDistance(dflag), alpha->symbol[i]);
    }
  }

  printf("\n");
  
  /* ---------------------------- */    
  /* now go ahead !               */

  skew = 0.;  /* just for SELECTOR_CUMUL */
  
  for (pos = 0 ; pos < len ; pos += step) {

    /* compute skew */
  
    if ((pos == 0) || (window < step)) {
      countingFunction(seq, pos, pos + window, allRight, nbins, 0);
      countingFunction(seq, pos + len - window, pos + len, allLeft, nbins, 0);
    }
    else {
      countingFunction(seq, pos - step, pos, allRight, nbins, -1);
      countingFunction(seq, pos - step + window, pos + window, allRight, nbins, 1);
      countingFunction(seq, pos - step + len - window, pos + len - window, allLeft, nbins, -1);
      countingFunction(seq, pos - step, pos, allLeft, nbins, 1);
    }
   
    mixingFunction(allLeft, allRight, countLeft, countRight, nbins);
   

    /* printout result */
    
    printf("%d", pos+1);

    /* printout skew */

    if (selector == SELECTOR_COUNT) {
    
      /* -- printout counts -- */
      
      for (i = 0 ; i <= alpha->length ; i++) {
        printf(" %d %d", countLeft[i], countRight[i]);
      }

      /* -- printout total alphabet chisquare -- */
    
      skew = distanceFunction(countLeft, countRight, nbins, -1);
      sign = ((nbins == 2) ? SignTwoBins(countRight, countLeft) : 1);
      printf(" %g", skew * sign);

      /* -- printout anti alphabet chisquare -- */

      tmpLeft[1] = countLeft[alpha->length];
      tmpRight[1] = countRight[alpha->length];
      tmpLeft[0] = tmpRight[0] = 0;
      for (i = 0 ; i < alpha->length ; i++) {
        tmpLeft[0] += countLeft[i];
        tmpRight[0] += countRight[i];
      }
      skew = distanceFunction(tmpLeft, tmpRight, 2, -1);
      sign = SignTwoBins(tmpRight, tmpLeft);
      printf(" %g", skew * sign);

      /* -- printout partial chisquares -- */

      tmpLeft[0] = tmpRight[0] = 0;
      tmpLeft[1] = tmpRight[1] = 0;
      for (j = 0 ; j < alpha->length ; j++) {
        tmpLeft[1] += countLeft[j];
        tmpRight[1] += countRight[j];
      }

      for (i = 0 ; i < alpha->length ; i++) {
        skew = distanceFunction(countLeft, countRight, nbins, i);
        tmpLeft[0]   = countLeft[i];
        tmpRight[0]  = countRight[i];
        tmpLeft[1]  -= countLeft[i];
        tmpRight[1] -= countRight[i];
        sign = SignTwoBins(tmpRight, tmpLeft);
        tmpLeft[1]  += countLeft[i];
        tmpRight[1] += countRight[i];
        printf(" %g", skew * sign);
      }
    }

    if (selector == SELECTOR_CUMUL) {

      /* -- printout XY cumulatives -- */
      
      sum = (float) ((countRight[0] + countLeft[0]) + (countRight[1] + countLeft[1]));
      if (sum != 0.) {
        skew += (float) ((countRight[0] + countLeft[0]) - (countRight[1] + countLeft[1]))
                / (sum * (float) len) * CUMUL_NORM;
      }
      printf(" %g", skew);
    }

    if (selector == SELECTOR_GENE) {

      /* -- printout direct -- */

      printf(" %d %d", countLeft[0], countRight[0]);

      /* -- printout reverse -- */

      printf(" %d %d", countLeft[1], countRight[1]);

      /* -- printout chi direct -- */

      skew = distanceFunction(countLeft, countRight, nbins, -1);
      sign = SignTwoBins(countRight, countLeft);
      printf(" %g", skew * sign);
    }

    if (selector == SELECTOR_CODON) {
    
      /* -- printout total codons -- */

      printf(" %d %d", countLeft[3* NB_NUC], countRight[3 * NB_NUC]);

      /* -- printout nuc-enrichment -- */

      for (i = 0 ; i < NB_NUC ; i++) {
        for (j = 0 ; j < 3 ; j++) {
          printf(" %d", countLeft[3*i + j]);
        }
        for (j = 0 ; j < 3 ; j++) {
          printf(" %d", countRight[3*i + j]);
        }
      }

      /* -- printout total chi nuc-enrichment -- */

      for (i = 0 ; i < NB_NUC ; i++) {
        skew = distanceFunction(countLeft + 3*i, countRight + 3*i, 3, -1);
        printf(" %g", skew);
      }
    }

    printf("\n");
  }

  FreeAlphabet(alpha);
  FreeSequence(seq);
  FreeZones(zones);

  exit(0);

}


        


