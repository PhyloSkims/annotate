/* ---------------------------------------------------------------- */
/* Copyright (c) Atelier de BioInformatique			    */
/* @file: util_complinv.c					    */
/* @desc: complinv a fasta sequence				    */
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
/*ARGSUSED*/

main(argn, argv)
	int  argn;
	char *argv[];
{
	int 	      nbseq, opt;
	FastaSequence *seq;

        while ((opt = getopt(argn, argv, "h")) != -1) {

	    switch (opt) {

		case 'h':
		    (void) printf("complement invert fasta sequence[s]\n");
		    (void) printf("usage: complinv\n");
		    exit(0);
		    break;

		case '?':
		    (void) printf("usage: complinv\n");
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
		continue;
	    }

	    (void) bio_seq_complement(seq->seq);
	    (void) str_reverse_string(seq->seq);
	    
	    (void) strcat(seq->comment, " (complement)");
	    
	    WriteFastaSequence(stdout, seq, FASTA_CHAR_PER_LINE);

	}

	FreeFastaSequence(seq);

	exit(0);
	
	/*NOTREACHED*/
}


	


