/* ---------------------------------------------------------------- */
/* Copyright (c) Atelier de BioInformatique			    */
/* @file: util_tab_nuc.c					    */
/* @desc: tabulate nucleotide usage				    */
/*								    */
/* @history:							    */
/* @+       <Gloup> : Jan 96 : PWG version			    */
/* ---------------------------------------------------------------- */

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef SGI
#include <getopt.h>
#endif

#include "Genetic.h"
#include "libbio.h"
#include "libfasta.h"

#ifndef Max
#define Max(i, j)  ((i) > (j) ? (i) : (j))
#endif

#define NB_NUC	 4

static char sDna[] = "ACTG";

/* ----------------------------------------------- */
static int sCharIndex(char *s, int c)
{
    char *ss = strchr(s, c);
    return (ss ? (int) (ss - s) : -1);
}

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

main(argn, argv)
	int  argn;
	char *argv[];
{
	int	      i, j, k, imin, imax, kmax, nbseq, opt;
	int	      sa_flag, so_flag, r_flag, t_flag, p_flag;
	int	      count[3][NB_NUC], tot[3][NB_NUC], sum[3];
	FastaSequence *seq;

	sa_flag = 0;		/* consider first codon		*/
	so_flag = 0;		/* consider last codon		*/
	t_flag  = 0;		/* no total			*/
	r_flag  = 0;		/* compute counts		*/
	p_flag  = 0;		/* no position specific		*/

        while ((opt = getopt(argn, argv, "hrsStp")) != -1) {
	    
	    switch (opt) {

		case 'h':
		    (void) printf("tabulate nucleotide usage\n");
		    (void) printf("usage: tab_nuc [-h] [-r] [-s] [-S] [-t] [-p]\n");
		    (void) printf("   -p       : counts by positions\n");
		    (void) printf("   -r       : compute relative frequencies\n");
		    (void) printf("   -s       : ignore first (start) codon\n");
		    (void) printf("   -S       : ignore last  (stop)  codon\n");
		    (void) printf("   -t       : print last total line\n");
		    exit(0);
		    break;

		case 't':
		    t_flag = 1;
		    break;

		case 'r':
		    r_flag = 1;
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
		    (void) printf("usage: tab_codon tab_nuc [-h] [-r] [-s] [-S] [-t] [-p]\n");
		    exit(6);
		    break;
	     }
	}

	seq = NewFastaSequence();

	nbseq = 0;

	kmax = (p_flag ? 3 : 1);
		    
	for (k = 0 ; k < kmax ; k++)
	for (i = 0 ; i < NB_NUC ; i++)
	    tot[k][i] = 0;

	printf("name");
	for (k = 0 ; k < kmax ; k++)
	for (i = 0 ; i < NB_NUC ; i++) {
	    if (p_flag)
		printf("\t%c%1d", sDna[i], (k+1));
	    else
		printf("\t%c", sDna[i]);
	}
	printf("\n");


	while (ReadFastaSequence(stdin, seq)) {

            nbseq++;

	    if (! seq->ok)
		(void) printf("error at seq # %d\n", nbseq);

	    for (k = 0 ; k < kmax ; k++)
	    for (i = 0 ; i < NB_NUC ; i++)
		count[k][i] = 0;
	    
	    imin  = (sa_flag ? 3 : 0);
	    imax  = seq->length - (so_flag ? 3 : 0);

	    for (i = imin, k = 0 ; i < imax ; i++) {

		j = sCharIndex(sDna, seq->seq[i]);

		if (j >= 0)
		    count[k][j]++;
		else
		    fprintf(stderr, "invalid nucleotide %1.1s at position %d in sequence %s\n",
				    seq->seq + i, i+1, seq->name);
		
		k = (p_flag ? (k+1)%3 : 0);
	    }

	    for (k = 0 ; k < kmax ; k++)
		sum[k] = 0;
	    if (r_flag) {
		for (k = 0 ; k < kmax ; k++)
		for (i = 0 ; i < NB_NUC ; i++)
		    sum[k] += count[k][i];
		for (k = 0 ; k < kmax ; k++)
		    sum[k] = Max(sum[k], 1);
	    }

	    printf("%s", seq->name);
	    
	    for (k = 0 ; k < kmax ; k++)
	    for (i = 0 ; i < NB_NUC ; i++) {
		if (r_flag)
		    printf("\t%.2f", (float) count[k][i] / (float) sum[k]);
		else
		    printf("\t%d", count[k][i]);
	    }
	    printf("\n");

	    if (t_flag) {
		for (k = 0 ; k < kmax ; k++)
		for (i = 0 ; i < NB_NUC ; i++)
		    tot[k][i] += count[k][i];
	    }

	}

	if (t_flag) {
	    printf("total:");

	    for (k = 0 ; k < kmax ; k++)
		sum[k] = 0;
	    if (r_flag) {
		for (k = 0 ; k < kmax ; k++)
		for (i = 0 ; i < NB_NUC ; i++)
		    sum[k] += tot[k][i];
		for (k = 0 ; k < kmax ; k++)
		    sum[k] = Max(sum[k], 1);
	    }

	    for (k = 0 ; k < kmax ; k++)
	    for (i = 0 ; i < NB_NUC ; i++) {
		if (r_flag)
		    printf("\t%.2f", (float) tot[k][i] / (float) sum[k]);
		else
		    printf("\t%d", tot[k][i]);
	    }
	    printf("\n");

	}	    

	FreeFastaSequence(seq);

	exit(0);
	
	/*NOTREACHED*/
}


	


