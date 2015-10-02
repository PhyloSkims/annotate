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
#define X_CODON  (NB_CODONS-1)
#define NB_AA	  21
#define X_AA	 (NB_AA-1)

typedef char NameString[16];

static char sDna[] = "ACTG";
static char sAA[]  = "ACDEFGHIKLMNPQRSTVWY*";

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

static void sMakeStdCodons(char *buffer)
{
    int i, j, k;
    
    for (i = 0 ; i < 4 ; i++)
    for (j = 0 ; j < 4 ; j++)
    for (k = 0 ; k < 4 ; k++) {
	*buffer++ = sDna[i];
	*buffer++ = sDna[j];
	*buffer++ = sDna[k];
	*buffer++ = '/';
    }
    *buffer = '\000';  
}

/* ----------------------------------------------- */
/* make string buffer of synonymous codons coding  */
/* for aa                                          */

static void sMakeAaCodons(char *buffer, int aa, int code)
{
    int	    aai;
    char    *c, codons[4*NB_CODONS+1];

    sMakeStdCodons(codons);  

    for (c = codons ; *c ; c += 4) {
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
/* return index of codon (range [0, X_CODONS])     */
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
/* or -1 if not found                     */

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
/* skew scores of codon                            */

static float sCodonScore(char *codon, int symb, int code) {
  int iaa, nsyn, n0, nM;
  char synon[4*NB_CODONS+1], buf[4];
  
  strncpy(buf, codon, 3);
  buf[3] = '\000';

  iaa = sAaIndex(buf, code);

  sMakeAaCodons(synon, sAA[iaa], code);
  
  n0 = sCount(buf, symb);  
  
  nsyn = strlen(synon) / 4;

  nM = sCount(synon, symb);

  return ((float) (n0) - ((float) (nM) / (float) nsyn));
}

/* ----------------------------------------------- */

main(argn, argv)
	int  argn;
	char *argv[];
{
	int	      i, j, k, imin, imax, nbseq, code, ncod;
	int	      opt, sa_flag, so_flag, p_flag;
	float         sscore[4], score[4][NB_CODONS];
	FastaSequence *seq;
	char	      *trip;
	char	      codons[4*NB_CODONS+1];

	extern char *optarg;	/* externs for getopts (3C)	*/

	code    = 0;		/* universal genetic code	*/
	sa_flag = 0;		/* consider first codon		*/
	so_flag = 0;		/* consider last codon		*/
	p_flag  = 0;		/* no pretty print		*/
	
	sMakeStdCodons(codons);	

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

	/* ------------------------------- */
	/* precompute score for each codon */
	/* ------------------------------- */

        for (j = 0 ; j < 4 ; j++) {
          for (i = 0, trip = codons ; i < X_CODON ; i++, trip += 4) {
            score[j][i] = sCodonScore(trip, sDna[j], code);
          }
        }

	seq = NewFastaSequence();

	nbseq = 0;

	if (p_flag) {
	    for (trip = codons ; *trip ; trip += 4) {
		printf("%3.3s", trip);
		k = sAaIndex(trip, code);
		printf("/%1c", sAA[k]);
		k = sCodonIndex(trip);
		for (j = 0; j < 4 ; j++)
		  printf(" %6.2f", score[j][k]);
		printf("\n");
	    }
	    printf("\n");
	}

	/* ---------------------------- */
	/* loop on sequences		*/
	/* ---------------------------- */

	while (ReadFastaSequence(stdin, seq)) {

            nbseq++;

	    if (! seq->ok)
		(void) printf("error at seq # %d\n", nbseq);

	    /* -------------------------------- */
	    /* compute score			*/
	    /* -------------------------------- */

	    for (j = 0 ; j < 4 ; j++)
              sscore[j] = 0.;

            ncod = 0;
            
	    imin  = (sa_flag ? 3 : 0);
	    imax  = seq->length - (so_flag ? 3 : 0);

	    for (i = imin ; i < imax ; i += 3) {

		k = sCodonIndex(seq->seq + i);

		if (k >= 0) {
		    for (j = 0 ; j < 4 ; j++) {
		      sscore[j] += (score[j][k] > 0 ? 1 : (score[j][k] < 0 ? -1 : 0));
		    }
		    ncod++;
		}
		else
		    fprintf(stderr, "invalid codon %3.3s at position %d in sequence %s\n",
				    seq->seq + i, i+1, seq->name);
	    }
	    
	    for (j = 0 ; j < 4 ; j++)
	      sscore[j] /= (float) ncod;
	      
            seq->comment[30] = '\000';
	    for (j = 0 ; j < 4 ; j++)
	      printf("%6.3f ", sscore[j]);
	    printf(" %s %s\n", seq->name, seq->comment);

	}
	
	/* ---------------------------- */
	/* end of read loop		*/
	/* ---------------------------- */

	/* ---------------------------- */
	/* free memory			*/
	/* ---------------------------- */

	FreeFastaSequence(seq);

	exit(0);
	
	/*NOTREACHED*/
}


	


