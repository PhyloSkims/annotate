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

#define COUNT_MODE  0	    /* comptage absolu	    	*/
#define FREQU_MODE  1	    /* freq. relative	    	*/
#define SYNON_MODE  2	    /* freq. rel. sur synonymes	*/

#define MIN_STORAGE 100

#define USE_HASH    1	    /* quicker version 		*/
			    /* ! see listing !		*/

typedef char NameString[16];

typedef struct s_Storage {
    NameString  name;
    int 	counts[NB_CODONS];
} Storage;


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
static char *sUpperProt(char *s)
{
    char *c;

    for (c = sUpper(s) ; *c ; c++) {
	if (*c == '#')
	    *c = '*';	
    }

    return s;
}


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
}

/* ----------------------------------------------- */
static void sMakeAaCodons(char *buffer, char *aa, char *codons, int code)
{
    int	    i, k;
    char    *c;

    while (*aa) {
	for (i = 0, c = codons ; i < X_CODON ; i++, c += 4) {
	    k = bio_codon_translate(c, code);
	    if (strchr(aa, k)) {
		strncpy(buffer, c, 3);
		buffer[3] = '/';
		buffer += 4;
	    }
	}
	aa++;
    } 
     
    *buffer = '\000';
}

/* ----------------------------------------------- */
static int sCodonIndex(char *s, char *c)
{
    int	i, ssa;
    char *sa, *sb;
    
    for (i = 0 ; i < 3 ; i++)
	if (! c[i]) return -1;
	
    sa  = c + 3;
    ssa = *sa;    
    *sa = '\000';
    sb = strstr(s, c);
    *sa = ssa;

    return (sb ? (int) (sb - s) / 4 : -1);
}

#if USE_HASH

/* ----------------------------------------------- */
/* this is a quicker version of sCodonIndex	   */
/* ! only works if alpha == sDna		   */

static int sCodonHash(char *alpha, char *c)
{
    int	i, h;
    char *p;
    
    for (i = h = 0 ; i < 3 ; i++) {
	if (! (p = strchr(alpha, c[i])))
	   return -1;
	h  = (h << 2) | (int) (p - alpha);
    }

    return h;
}

#endif

/* ----------------------------------------------- */
static int sAaIndex(char *saa, char *c, int code)
{
    int  aa;
    char *paa;
    
    aa = bio_codon_translate(c, code);
    if ((paa = strchr(saa, aa)) != 0) 
	return (int) (paa - saa);
    return X_AA;
}

/* ----------------------------------------------- */
static int sSumCounts(char *codons, int *counts, int zero)
{
    int  i, sum;
    char *trip;

    for (i = sum = 0, trip = codons ; i < X_CODON ; i++, trip += 4)
	sum += counts[i];

    return (sum ? sum : zero);
}

/* ----------------------------------------------- */
static void sSumAA(char *codons, char *aa, int code, int *counts, int *naa, int zero)
{
    int  i, k;
    char *trip;

    for (i = 0 ; i < NB_AA ; i++)
	naa[i] = 0;

    for (i = 0, trip = codons ; i < X_CODON ; i++, trip += 4) {
	k = sAaIndex(aa, trip, code);
	naa[k] += counts[i];
    }
    
    for (i = 0 ; i < NB_AA ; i++)
	if (naa[i] == 0)
	   naa[i] = zero;
}   


/* ----------------------------------------------- */
static void sPrintCodonUsage(int *tot, char *codons, int code)
{
    int	    k, n;
    char    *aa, *trip;
    char    bufaa[2],
	    bufco[4*NB_CODONS+1];

    bufaa[1] = '\000';
    
    for (aa = sAA ; *aa ; aa++) {
	*bufaa = *aa;
	sMakeAaCodons(bufco, bufaa, codons, code);
	printf ("%1s\n", bufaa);
	for (trip = bufco, n = 0 ; *trip ; trip += 4) {
	    k = sCodonIndex(codons, trip);
	    n += tot[k];
	}
	if (n == 0) n = 1;
	for (trip = bufco ; *trip ; trip += 4) {
	    k = sCodonIndex(codons, trip);
	    printf("  %3.3s\t%.2f %d\n", trip, (float) tot[k]/ (float) n,
					   tot[k]);
	}
    }
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
	int	      i, j, k, imin, imax, nbseq, opt, code, sum, npart, total;
	int	      sa_flag, so_flag, r_flag, t_flag, p_flag, g_flag, a_flag;
	int	      counts[NB_CODONS], totcd[NB_CODONS], corcd[NB_CODONS],
		      naa[NB_AA], totaa[NB_AA];
	Storage	      *partial;	      
	FastaSequence *seq;
	char	      *trip;
	char	      codons[4*NB_CODONS+1],
		      aas_ignore[256], 
		      usr_ignore[4*NB_CODONS+1]; 

	extern char *optarg;	/* externs for getopts (3C)	*/

	code    = 0;		/* universal genetic code	*/
	sa_flag = 0;		/* consider first codon		*/
	so_flag = 0;		/* consider last codon		*/
	t_flag  = 0;		/* no total			*/
	r_flag  = COUNT_MODE;	/* compute counts		*/
	p_flag  = 0;		/* no pretty print		*/
	g_flag  = 0;		/* no global correction		*/
	a_flag  = 0;		/* no aa names			*/
	
	*aas_ignore = '\000';
	*usr_ignore = '\000';
	
	npart   = 0;
	partial = NULL;
	
	sMakeStdCodons(codons);	

	/* ---------------------------- */
	/* parse arguments 		*/
	/* ---------------------------- */
			
        while ((opt = getopt(argn, argv, "ac:hgi:I:rRsStp")) != -1) {
	    
	    switch (opt) {

		case 'a':
		    a_flag = 1;
		    break;

		case 'c':
		    if (   (sscanf(optarg, "%d", &code) != 1)
			|| (code < 0) || (code > 8)) {
		       (void) printf("bad code value: -c (0-8)\n");
		       exit(5);
		    }
		    break;

		case 'h':
		    (void) printf("tabulate codon usage (CU)\n");
		    (void) printf("usage: tab_codon [-a] [-c code] [-i|I alpha]\n");
		    (void) printf("                 [-r|R] [-s] [-S] [-t] [-p]\n");
		    (void) printf("   -a   : add the amino-acid symbol in axis name\n");
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
		    (void) printf("   -i alpha : ignore codons in alpha\n");
		    (void) printf("      alpha has form: 'ATG/GTG/CTG'\n");
		    (void) printf("   -I alpha : ignore codons whose aa are in alpha\n");
		    (void) printf("      alpha has form: 'FYW'\n");
		    (void) printf("                      '#' means stop\n");
		    (void) printf("   -r       : compute relative frequencies (rCU)\n");
		    (void) printf("   -R       : compute synonymous relative frequencies (rSCU)\n");
		    (void) printf("              note: by default (no -r nor -R) counts (CU) are printed\n");
		    (void) printf("   -g       : apply various global count corrections :\n");
		    (void) printf("              for CU   : global aa usage correction -> SCU\n");
		    (void) printf("              for rCU  : global aa usage correction -> SrCU\n");
		    (void) printf("              for rSCU : global codon usage correction ->CrSCU\n");
		    (void) printf("   -s       : ignore first (start) codon\n");
		    (void) printf("   -S       : ignore last  (stop)  codon\n");
		    (void) printf("   -p       : pretty print codon usage\n");
		    (void) printf("              -r -R -g options are then ignored\n");
		    (void) printf("   -t       : print last total line\n");
		    exit(0);
		    break;

		case 'g':
		    g_flag = 1;
		    break;

		case 'i':
		    (void) strcpy(usr_ignore, optarg);
		    (void) sUpperDna(usr_ignore);
		    break;

		case 'I':
		    (void) strcpy(aas_ignore, optarg);
		    (void) sUpperProt(aas_ignore);
		    break;

		case 't':
		    t_flag = 1;
		    break;

		case 'r':
		    r_flag = FREQU_MODE;
		    break;

		case 'R':
		    r_flag = SYNON_MODE;
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
		    (void) printf("usage: tab_codon [-h] [-c code] [-i|I alpha]\n");
		    (void) printf("                 [-r|R] [-s] [-S] [-t] [-p] [-a] [-b]\n");
		    exit(6);
		    break;
	     }
	}

	/* ---------------------------- */
	/* check usage			*/
	/* ---------------------------- */

	if (*aas_ignore && *usr_ignore) {
		fprintf(stderr,"tab_codon: -i and -I incompatible options\n");
		exit(5);
	}

	if (*aas_ignore)
	    sMakeAaCodons(usr_ignore, aas_ignore, codons, code);

	seq = NewFastaSequence();

	nbseq = 0;

	for (i = 0 ; i < NB_CODONS ; i++)
	    totcd[i] = 0;

	if (! p_flag) {
	    printf("name");
	    for (i = 0, trip = codons ; i < X_CODON ; i++, trip += 4) {
		if (sCodonIndex(usr_ignore, trip) == -1) {
		    printf("\t%3.3s", trip);
		    if (a_flag) {
		       k = sAaIndex(sAA, trip, code);
		       printf("/%1c", sAA[k]);
		    }
		}
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
	    /* compute counts			*/
	    /* -------------------------------- */

	    for (i = 0 ; i < NB_CODONS ; i++)
		counts[i] = 0;
	    
	    imin  = (sa_flag ? 3 : 0);
	    imax  = seq->length - (so_flag ? 3 : 0);

  	    k = -1;
	    
	    for (i = imin ; i < imax ; i += 3) {

#if USE_HASH
		k = sCodonHash(sDna, seq->seq + i);
#else
		k = sCodonIndex(codons, seq->seq + i);
#endif

		if (k >= 0)
		    counts[k]++;
		else
		    fprintf(stderr, "invalid codon %3.3s at position %d in sequence %s\n",
				    seq->seq + i, i+1, seq->name);
	    }

	    /* remove ignored codons		*/

	    for (i = 0, trip = codons ; i < X_CODON ; i++, trip += 4) {

		if (sCodonIndex(usr_ignore, trip) >= 0)
		    counts[i] = 0;
	    }

	    /* -------------------------------- */
	    /* compute totals			*/
	    /* -------------------------------- */
	    
	    for (i = 0 ; i < NB_CODONS ; i++)
		totcd[i] += counts[i];

	    /* -------------------------------- */
	    /* store or print local values	*/
	    /* -------------------------------- */

	    if (! p_flag) {
					/* ------------ */
		if (g_flag) {		/* store	*/
					/* ------------ */

	    	    if (npart <= nbseq) {
		       partial = sIncreaseStorage(partial, &npart);
		       if (! partial) {
			    fprintf(stderr,"not enough memory for %d sequences\n", nbseq);
			    exit(10);
		       }
		    }

		    sCopyStorage(partial + nbseq - 1, counts, seq->name);

		}

					/* ------------ */
		else {			/* printout	*/
					/* ------------ */
		    sum = sSumCounts(codons, counts, 1);

		    if (r_flag == SYNON_MODE)
			sSumAA(codons, sAA, code, counts, naa, 1);

		    printf("%s", seq->name);

		    for (i = 0, trip = codons ; i < X_CODON ; i++, trip += 4) {
			if (sCodonIndex(usr_ignore, trip) == -1) {
			    if (r_flag == COUNT_MODE)
				printf("\t%d", counts[i]);
			    else if (r_flag == FREQU_MODE)
				printf("\t%.2f", 100. * (float) counts[i] / (float) sum);
			    else if (r_flag == SYNON_MODE) {
				k = sAaIndex(sAA, trip, code);
				printf("\t%.1f", 100. * (float) counts[i] / (float) naa[k]);
			    }
			}
		    }
		    printf("\n");
		}
	    }

	}
	
	/* ---------------------------- */
	/* end of read loop		*/
	/* ---------------------------- */

	/* ---------------------------- */
	/* global correction		*/
	/* ---------------------------- */

	if (g_flag && (! p_flag)) {

	    /* ------------------------ */
	    /* compute corrected counts	*/
	    /* ------------------------ */

	    total = sSumCounts(codons, totcd, 1);

	    sSumAA(codons, sAA, code, totcd, totaa, 1);

	    if ((r_flag == COUNT_MODE) || (r_flag == FREQU_MODE)) {

	    	sSumAA(codons, sAA, code, totcd, totaa, 0);

		for (i = 0 ; i < X_CODON ; i++)
		    corcd[i] = 0;
	    	for (j = 0 ; j < nbseq ; j++) {
		   sum = sSumCounts(codons, partial[j].counts, 0);
		   sSumAA(codons, sAA, code, partial[j].counts, naa, 1);
	    	   for (i = 0, trip = codons ; i < X_CODON ; i++, trip += 4) {
			k = sAaIndex(sAA, trip, code);
			partial[j].counts[i] = (int) floor(  
						      (double)partial[j].counts[i] * (double) totaa[k]
		       				    * (double) sum / (double) total / (double) naa[k]
						    + 0.5);
			corcd[i] += partial[j].counts[i];
		   }
		}
	    }
	    else {
	    	for (i = 0 ; i < X_CODON ; i++)
		    corcd[i] = totcd[i];
	    
	    }
	    	
	    
	    /* ------------------------ */
	    /* printout sequences	*/
	    /* ------------------------ */

	    
	    for (j = 0 ; j < nbseq ; j++) {

		sum = sSumCounts(codons, partial[j].counts, 1);

		if (r_flag == SYNON_MODE)
		    sSumAA(codons, sAA, code, partial[j].counts, naa, 0);
		
		printf("%s", partial[j].name);

		for (i = 0, trip = codons ; i < X_CODON ; i++, trip += 4) {
		    if (sCodonIndex(usr_ignore, trip) == -1) {
			if (r_flag == COUNT_MODE)
			    printf("\t%d", partial[j].counts[i]);
			else if (r_flag == FREQU_MODE)
			    printf("\t%.2f", 100. * (float) partial[j].counts[i] / (float) sum);
			else if (r_flag == SYNON_MODE) {
			    k = sAaIndex(sAA, trip, code);
			    if (naa[k] != 0)
			    	printf("\t%.1f", 100. * (float) partial[j].counts[i] / (float) naa[k]);
			    else
			    	printf("\t%.1f", 100. * (float) totcd[i] / (float) totaa[k]);
			}
		    }
		}
		printf("\n");
	    }

	    total = sSumCounts(codons, corcd, 1);
	    
	    for (i = 0 ; i < X_CODON ; i++)
		totcd[i] = corcd[i];

	    sSumAA(codons, sAA, code, totcd, totaa, 1);

	}

	/* ---------------------------- */
	/* print grand total		*/
	/* ---------------------------- */

	if (t_flag && (! p_flag)) {

	    printf("total:");

	    sum = sSumCounts(codons, totcd, 1);
	    
	    for (i = 0, trip = codons ; i < X_CODON ; i++, trip += 4) {
		if (sCodonIndex(usr_ignore, trip) == -1) {
		    if (r_flag == COUNT_MODE)
			printf("\t%d", totcd[i]);
		    else if (r_flag == FREQU_MODE)
			printf("\t%.2f", 100. * (float) totcd[i] / (float) sum);
		    else if (r_flag == SYNON_MODE) {
			k = sAaIndex(sAA, trip, code);
			printf("\t%.1f", 100. * (float) totcd[i] / (float) totaa[k]);
		    }
		}
	    }
	    printf("\n");

	}

	/* ---------------------------- */
	/* pretty print			*/
	/* ---------------------------- */

	if (p_flag) {
	    sPrintCodonUsage(totcd, codons, code);
	}

	/* ---------------------------- */
	/* free memory			*/
	/* ---------------------------- */

#if 1

	if (partial)
	   free(partial);

	FreeFastaSequence(seq);

#endif

	exit(0);
	
	/*NOTREACHED*/
}


	


