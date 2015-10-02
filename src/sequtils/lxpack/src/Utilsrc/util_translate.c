/* ---------------------------------------------------------------- */
/* Copyright (c) Atelier de BioInformatique			    */
/* @file: util_translate.c					    */
/* @desc: util_translate a fasta nucleic sequence		    */
/*								    */
/* @history:							    */
/* @+       <Gloup> : Jan 96 : PWG version			    */
/* ---------------------------------------------------------------- */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef SGI
#include <getopt.h>
#endif

#include "libbio.h"
#include "libfasta.h"



/* ----------------------------------------------- */
/* safely copy overlapping strings		   */

static void sStrCpy(char *s1, char *s2)
{
    while (*s2)
       *s1++ = *s2++;
    *s1 = '\000';
}

/* ----------------------------------------------- */

main(argn, argv)
	int  argn;
	char *argv[];
{
	int	      nbseq, minlen, opt, code;
	int	      sa_flag, so_flag;
	FastaSequence *seq;

	extern char *optarg;	/* externs for getopts (3C)	*/

	code = 0;
	sa_flag = so_flag = 0;
	
        while ((opt = getopt(argn, argv, "c:hsS")) != -1) {
	    
	    switch (opt) {

		case 'c':
		    if (   (sscanf(optarg, "%d", &code) != 1)
			|| (code < 0) || (code > 8)) {
		       (void) printf("bad code value: -c (0-8)\n");
		       exit(5);
		    }
		    break;

		case 'h':
		    (void) printf("translate dna sequence[s] to protein[s]\n");
		    (void) printf("usage: translate [-c (0-8)] [-s] [-S]\n");
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
		    (void) printf("   -s : ignore first (start) codon\n");
		    (void) printf("   -S : ignore last  (stop)  codon\n");
		    exit(0);
		    break;

		case 's':
		    sa_flag = 1;
		    break;

		case 'S':
		    so_flag = 1;
		    break;

		case '?':
		    (void) printf("usage: translate [-c (0-8)] [-s] [-S]\n");
		    exit(6);
		    break;
	     }
	}

	seq = NewFastaSequence();

	nbseq = 0;

	while (ReadFastaSequence(stdin, seq)) {

            nbseq++;

	    if (! seq->ok) {
		(void) printf("error at seq # %d\n", nbseq);
		(void) printf("bad length at seq # %d\n", nbseq);
		continue;
	    }

	    minlen = 0;
	    if (sa_flag) minlen += 3;
	    if (so_flag) minlen += 3;

	    if (seq->length <= minlen) {
		(void) printf("bad length at seq # %d\n", nbseq);
		continue;
	    }

	    (void) bio_seq_translate(seq->seq, code);

	    (void) strcat(seq->comment, " (translation)");

	    seq->length /= 3;
	    
	    if (sa_flag) {
		sStrCpy(seq->seq, seq->seq + 1);
		seq->length--;
	    }
	    
	    if (so_flag)
		seq->length--;
	    
	    WriteFastaSequence(stdout, seq, FASTA_CHAR_PER_LINE);

	}

	FreeFastaSequence(seq);

	exit(0);
	
	/*NOTREACHED*/
}


	


