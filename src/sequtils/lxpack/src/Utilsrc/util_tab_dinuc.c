/* ---------------------------------------------------------------- */
/* Copyright (c) Atelier de BioInformatique			    */
/* @file: util_tab_nuc.c					    */
/* @desc: tabulate dinucleotide usage				    */
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

#define NB_DINUC 16

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
static int sDiHash(char *alpha, char *c)
{
    int	i, h;
    char *p;
    
    for (i = h = 0 ; i < 2 ; i++) {
	if (! (p = strchr(alpha, c[i])))
	   return -1;
	h  = (h << 2) | (int) (p - alpha);
    }

    return h;
}


/* ----------------------------------------------- */
static int sUpdateDiHash(int oldh, char *alpha, char *c)
{
    char *p;

    if (oldh < 0)
       return sDiHash(alpha, c);

    if (! (p = strchr(alpha, c[1])))
	return -1;
	
    return (((oldh << 2) & 0xf) |  (int) (p - alpha));
}

/* ----------------------------------------------- */
static void sMakeStdDinuc(char *buffer)
{
    int i, j;
    
    for (i = 0 ; i < 4 ; i++)
    for (j = 0 ; j < 4 ; j++) {
	*buffer++ = sDna[i];
	*buffer++ = sDna[j];
	*buffer++ = '/';
    }    
}

/* ----------------------------------------------- */
static int sInPhase(int phase, int flag)
{
    switch (flag) {
	case 1   : return  (phase == 0);
	case 2   : return  (phase == 1);
	case 3   : return  (phase == 2);
	case 12  : return ((phase == 0) || (phase == 1));
	case 13  : return ((phase == 0) || (phase == 2));
	case 23  : return ((phase == 1) || (phase == 2));
	case 123 : return 1;
   }

   return 0;
}

/* ----------------------------------------------- */

main(argn, argv)
	int  argn;
	char *argv[];
{
	int	      i, j, k, imin, imax, sum, nbseq, opt;
	int	      sa_flag, so_flag, r_flag, t_flag, p_flag;
	int	      count[NB_DINUC], tot[NB_DINUC];
	FastaSequence *seq;
	char	      *din, dinuc[3*NB_DINUC+1];
	
	extern char *optarg;	/* externs for getopts (3C)	*/

	sa_flag = 0;		/* consider first codon		*/
	so_flag = 0;		/* consider last codon		*/
	t_flag  = 0;		/* no total			*/
	r_flag  = 0;		/* compute counts		*/
	p_flag  = 123;		/* 3 phases			*/

	sMakeStdDinuc(dinuc);	

        while ((opt = getopt(argn, argv, "hrsStp:")) != -1) {
	    
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
		    if (    (sscanf(optarg, "%d", &p_flag) != 1)
		        || (   (p_flag != 1)  && (p_flag != 2)  && (p_flag != 3)
			    && (p_flag != 12) && (p_flag != 13) && (p_flag != 23)
			    && (p_flag != 123))) {
		       (void) printf("bad phase value: -p [1][2][3]\n");
		       exit(5);
		    }
		    break;
		   
		case '?':
		    (void) printf("usage: tab_codon tab_nuc [-h] [-r] [-s] [-S] [-t] [-p [123]]\n");
		    exit(6);
		    break;
	     }
	}

	seq = NewFastaSequence();

	nbseq = 0;

	for (i = 0 ; i < NB_DINUC ; i++)
	    tot[i] = 0;

	printf("name");
	for (i = 0, din = dinuc ; i < NB_DINUC ; i++, din += 3) {
		printf("\t%2.2s", din);
	}
	printf("\n");

	while (ReadFastaSequence(stdin, seq)) {

            nbseq++;

	    if (! seq->ok)
		(void) printf("error at seq # %d\n", nbseq);

	    for (i = 0 ; i < NB_DINUC ; i++)
		count[i] = 0;
	    
	    imin  = (sa_flag ? 3 : 0);
	    imax  = seq->length - (so_flag ? 3 : 0) - 1;

	    for (i = imin, j = -1, k = 0 ; i < imax ; i++) {

		j = sUpdateDiHash(j, sDna, seq->seq + i);

		if (sInPhase(k, p_flag)) {
		    if (j >= 0)
			count[j]++;
		    else
			fprintf(stderr, "invalid dinucleotide %2.2s at position %d in sequence %s\n",
					seq->seq + i, i+1, seq->name);
		}
				
		k = (k+1) % 3;
	    }

	    sum = 0;
	    
	    if (r_flag) {
		for (i = 0 ; i < NB_DINUC ; i++)
		    sum += count[i];
		sum = Max(sum, 1);
	    }

	    printf("%s", seq->name);
	    
	    for (i = 0 ; i < NB_DINUC ; i++) {
		if (r_flag)
		    printf("\t%.2f", (float) count[i] / (float) sum);
		else
		    printf("\t%d", count[i]);
	    }
	    printf("\n");

	    if (t_flag) {
		for (i = 0 ; i < NB_DINUC ; i++)
		    tot[i] += count[i];
	    }

	}

	if (t_flag) {
	    printf("total:");

	    sum = 0;
	    if (r_flag) {
		for (i = 0 ; i < NB_DINUC ; i++)
		    sum += tot[i];
		sum = Max(sum, 1);
	    }

	    for (i = 0 ; i < NB_DINUC ; i++) {
		if (r_flag)
		    printf("\t%.2f", (float) tot[i] / (float) sum);
		else
		    printf("\t%d", tot[i]);
	    }
	    printf("\n");

	}	    

	FreeFastaSequence(seq);

	exit(0);
	
	/*NOTREACHED*/
}


	


