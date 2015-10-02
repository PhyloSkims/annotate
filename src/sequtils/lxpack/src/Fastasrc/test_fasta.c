/* ---------------- */
/* @file: ctest.c   */
/* ---------------- */

#include <stdio.h>
#include <stdlib.h>

#include "libfasta.h"

/* ----------------------------------------------- */

main(argn, argv)
	int  argn;
	char *argv[];
{
	int	      nbseq;
	FastaSequence *seq;


	seq = NewFastaSequence();

	nbseq = 0;

	while (ReadFastaSequence(stdin, seq)) {

	if (! seq->ok)
	printf("error\n");

            nbseq++;

	    if (nbseq % 1000 == 999) {
		printf("\r%d", nbseq);
		fflush(stdout);
	    }


	    WriteFastaSequence(stdout, seq, FASTA_CHAR_PER_LINE);

	}

	FreeFastaSequence(seq);

	printf("total seq = %d\n", nbseq);

	exit(0);

}


	


