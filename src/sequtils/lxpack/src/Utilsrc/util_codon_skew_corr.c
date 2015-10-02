/* ---------------------------------------------------------------- */
/* Copyright (c) Atelier de BioInformatique			    */
/* @file: util_tab_codon.c					    */
/* @desc: tabulate codon usage					    */
/*								    */
/* @history:							    */
/* @+       <Gloup> : Jan 96 : PWG version			    */
/* ---------------------------------------------------------------- */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifdef SGI
#include <getopt.h>
#endif

#include "Genetic.h"
#include "libbio.h"
#include "libfasta.h"

#ifndef Max
#define Max(i, j)  ((i) > (j) ? (i) : (j))
#endif

#define NB_CODONS 65
#define NB_AA	  21

#define MIN_STORAGE 1000

typedef char NameString[16];

typedef struct s_Storage {
    NameString  name;
    int 	counts[NB_CODONS];
} Storage;

static char sDna[] = "ACTG";
static char sAA[]  = "ACDEFGHIKLMNPQRSTVWY*";
static char sCodon[4*NB_CODONS+1];


/* ----------------------------------------------- */
static char *sUpper(char *s)
{
    char *c;

    for (c = s ; *c ; c++) {
	if (islower(*c))
	    *c = toupper(*c);
    }

    return s;
}

/* ----------------------------------------------- */
static char *sUpperDna(char *s)
{
    char *c;

    for (c = sUpper(s) ; *c ; c++) {
	if (*c == 'U')
	    *c = 'T';	
    }

    return s;
}

/* ----------------------------------------------- */
/* make string buffer of all codons :              */
/* "AAA/AAC/AAT/AAG/ACA/....../GGG/"               */
/* ----------------------------------------------- */

static void sMakeStdCodons()
{
    int i, j, k;
    char *buf = sCodon;
    
    for (i = 0 ; i < 4 ; i++)
    for (j = 0 ; j < 4 ; j++)
    for (k = 0 ; k < 4 ; k++) {
	*buf++ = sDna[i];
	*buf++ = sDna[j];
	*buf++ = sDna[k];
	*buf++ = '/';
    }
    *buf = '\000';  
}

/* ----------------------------------------------- */
/* make string buffer of synonymous codons coding  */
/* for aa                                          */

static void sMakeAaCodons(char *buffer, int aa, int code)
{
    int	    aai;
    char    *c;

    for (c = sCodon ; *c ; c += 4) {
      aai = bio_codon_translate(c, code);
      if (aai == aa) {
        strncpy(buffer, c, 3);
	buffer[3] = '/';
	buffer += 4;
      }
    }
    *buffer = '\000';
}

/* ----------------------------------------------- */
/* return index of codon (range [0, NB_CODONS[)    */
/* return -1 if a symbol is not found              */

static int sCodonIndex(char *codon)
{
    int	i, h;
    char *p;
    
    for (i = h = 0 ; i < 3 ; i++) {
	if (! (p = strchr(sDna, codon[i])))
	   return -1;
	h  = (h << 2) | (int) (p - sDna);
    }

    return h;
}

/* ----------------------------------------------- */
/* return index of aa encoded by codon             */
/* or -1 if not found                              */

static int sAaIndex(char *codon, int code)
{
    int  aa;
    char *paa;
    
    aa = bio_codon_translate(codon, code);
    if ((paa = strchr(sAA, aa)) != 0) 
	return (int) (paa - sAA);
    return -1;
}

/* ----------------------------------------------- */
/* count # occurences of char mark in buffer       */

static int sCount(char *buffer, int mark)
{
    int count;
    
    for (count = 0 ; *buffer ; buffer++) {
      if (*buffer == mark)
        count++;
    }
    
    return count;
}

/* ----------------------------------------------- */
/* compute relative synonymous frequency of codons */

static void sSynFrequency(int *count, float *freq, int code) {
  int icod, k, iaa, sum;
  char synon[4*NB_CODONS+1], buf[4];
  char *trip;

  buf[3] = '\000';  
  for (icod = 0 ; icod < (NB_CODONS-1) ; icod++) {
    strncpy(buf, sCodon + 4 * icod, 3);
    iaa = sAaIndex(buf, code);
    sMakeAaCodons(synon, sAA[iaa], code);
    sum = 0;
    for (trip = synon ; *trip ; trip += 4) {
      k = sCodonIndex(trip);
      sum += count[k];
    }
    freq[icod] = (float) count[icod] / (float) sum;
  }
}


/* ----------------------------------------------- */
/* skew score of codon                             */

static float sCodonScore(char *codon, float *freq, int code) {
  int iaa, icod;
  float xs0, xsM;
  char synon[4*NB_CODONS+1], buf[4];
  char *trip;
  
  strncpy(buf, codon, 3);
  buf[3] = '\000';

  iaa = sAaIndex(buf, code);

  sMakeAaCodons(synon, sAA[iaa], code);
  
  xs0 = (float) (sCount(buf, 'G') - sCount(buf, 'C'));
  
  xsM = 0.;
  for (trip = synon ; *trip ; trip += 4) {
    strncpy(buf, trip, 3);
    icod = sCodonIndex(buf);
    xsM += (float) (sCount(buf, 'G') - sCount(buf, 'C')) * freq[icod];
  }
  
  return (xs0 - xsM);
}

/* ----------------------------------------------- */
static Storage *sIncreaseStorage(Storage *store, int *size)
{
    int      nsiz;
    Storage *new;

    nsiz = Max(*size * 2, MIN_STORAGE);
    

    if (store)
    	new = (Storage *) realloc(store, nsiz * sizeof(Storage));
    else
        new = (Storage *) malloc(nsiz * sizeof(Storage));

    if (new)
       *size = nsiz;
       
    return new;
}

/* ----------------------------------------------- */
static void sCopyStorage(Storage *store, int *counts, char *name)
{
    int i;
    
    for (i = 0 ; i < NB_CODONS ; i++)
        store->counts[i] = counts[i];
	
    (void) strncpy(store->name, name, sizeof(NameString));
    store->name[sizeof(NameString)-1] = '\000';
}    

/* ----------------------------------------------- */

main(argn, argv)
	int  argn;
	char *argv[];
{
	int	      i, k, imin, imax, nbseq, code, nstore, sign, count, sum;
	int	      opt, sa_flag, so_flag, p_flag;
	int	      counts[NB_CODONS], totcd[NB_CODONS];
	float         freq[NB_CODONS], score[NB_CODONS], sscore;
	Storage	      *store;	      
	FastaSequence *seq;
	char	      *trip;

	extern char *optarg;	/* externs for getopts (3C)	*/

	code    = 0;		/* universal genetic code	*/
	sa_flag = 0;		/* consider first codon		*/
	so_flag = 0;		/* consider last codon		*/
	p_flag  = 0;		/* no pretty print		*/

	nstore   = 0;
	store    = NULL;
	
	sMakeStdCodons();	

	/* ---------------------------- */
	/* parse arguments 		*/
	/* ---------------------------- */
			
        while ((opt = getopt(argn, argv, "c:hsSp")) != -1) {
	    
	    switch (opt) {

		case 'c':
		    if (   (sscanf(optarg, "%d", &code) != 1)
			|| (code < 0) || (code > 8)) {
		       (void) printf("bad code value: -c (0-8)\n");
		       exit(5);
		    }
		    break;

		case 'h':
		    (void) printf("codon GC skew\n");
		    (void) printf("usage: codon_skew [-c code]\n");
		    (void) printf("                  [-s] [-S] [-p]\n");
		    (void) printf("   -c code\n");
		    (void) printf("      0 : universal\n");
		    (void) printf("      1 : mito yeast\n");
		    (void) printf("      2 : mito vertebrate\n");
		    (void) printf("      3 : filamentous fungi\n");
		    (void) printf("      4 : mito insects & platyhelminthes\n");
		    (void) printf("      5 : Candida cylindracea\n");
		    (void) printf("      6 : Ciliata\n");
		    (void) printf("      7 : Euplotes\n");
		    (void) printf("      8 : mito echinoderms\n");
		    (void) printf("   -s       : ignore first (start) codon\n");
		    (void) printf("   -S       : ignore last  (stop)  codon\n");
		    (void) printf("   -p       : pretty print codon score\n");
		    exit(0);
		    break;

		case 's':
		    sa_flag = 1;
		    break;

		case 'S':
		    so_flag = 1;
		    break;

		case 'p':
		    p_flag = 1;
		    break;
		   
		case '?':
		    (void) printf("usage: codon_skew [-h] [-c code]\n");
		    (void) printf("                  [-s] [-S] [-p]\n");
		    exit(6);
		    break;
	     }
	}

	/* ---------------------------- */
	/* check usage			*/
	/* ---------------------------- */

	/* ---------------------------- */
	/* loop on sequences		*/
	/* ---------------------------- */

	for (i = 0 ; i < NB_CODONS ; i++)
	  totcd[i] = 0;

	seq = NewFastaSequence();

	nbseq = 0;

	while (ReadFastaSequence(stdin, seq)) {

            nbseq++;

	    if (! seq->ok)
		(void) printf("error at seq # %d\n", nbseq);

	    /* -------------------------------- */
	    /* compute counts			*/
	    /* -------------------------------- */

	    for (i = 0 ; i < NB_CODONS ; i++)
		counts[i] = 0;

	    imin  = (sa_flag ? 3 : 0);
	    imax  = seq->length - (so_flag ? 3 : 0);

	    for (i = imin ; i < imax ; i += 3) {

		k = sCodonIndex(seq->seq + i);

		if (k >= 0)
		    counts[k]++;
		else
		    fprintf(stderr, "invalid codon %3.3s at position %d in sequence %s\n",
				    seq->seq + i, i+1, seq->name);
	    }

	    /* -------------------------------- */
	    /* compute grand total		*/
	    /* -------------------------------- */
	    
	    for (i = 0 ; i < NB_CODONS ; i++)
		totcd[i] += counts[i];
	    
	    /* -------------------------------- */
	    /* store counts			*/
	    /* -------------------------------- */
	    
	    if (nstore <= nbseq) {
	       store = sIncreaseStorage(store, &nstore);
	       if (! store) {
		    fprintf(stderr,"not enough memory for %d sequences\n", nbseq);
		    exit(10);
	       }
	    }

	    sCopyStorage(store + nbseq - 1, counts, seq->name);
	}

	/* ------------------------------- */
	/* compute score for each codon    */
	/* ------------------------------- */

	sSynFrequency(totcd, freq, code);
	
        for (i = 0, trip = sCodon ; *trip ; i++, trip += 4) {
          score[i] = sCodonScore(trip, freq, code);
        }

	if (p_flag) {
	    for (trip = sCodon ; *trip ; trip += 4) {
		printf("%3.3s", trip);
		k = sAaIndex(trip, code);
		printf("/%1c", sAA[k]);
		k = sCodonIndex(trip);
		printf(" %5.3f", freq[k]);
		printf(" %6.2f\n", score[k]);
	    }
	    printf("\n");
	}

	/* ------------------------------- */
	/* compute score for each sequence */
	/* ------------------------------- */

        for (i = 0 ; i < nbseq ; i++) {
          sscore = 0.;
          sum = 0;
          for (k = 0 ; k < NB_CODONS ; k++) {
            sign = (score[k] > 0 ? 1 : (score[k] < 0 ? -1 : 0));
            count = store[i].counts[k];
            sum += count;
            sscore += (float) sign * (float) count;
          }
          sscore /= (float) sum;
	  printf("%s %6.3f\n", store[i].name, sscore);
        }


	/* ---------------------------- */
	/* free memory			*/
	/* ---------------------------- */

	FreeFastaSequence(seq);
	
	free(store);

	exit(0);
	
	/*NOTREACHED*/
}


	


