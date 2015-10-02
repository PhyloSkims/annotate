/* ---------------------------------------------------------------- */
/* Copyright (c) Atelier de BioInformatique			    */
/* @file: util_cut.c						    */
/* @desc: cut a fasta sequence					    */
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
/* printout sequence				   */

static void sPrintSequence(FastaSequence *seq, int from, int to, int rev, char *name) {
    int ifrom, ito, ilen;
    
    static FastaSequence *bufSeq = NULL;

    if (bufSeq == NULL)
       bufSeq = NewFastaSequence();
    
    ifrom = ((from > 0) ? from : seq->length + from);
    
    ito = ((to > 0) ? to : seq->length + to);

    if ((ito > seq->length) || (ifrom > ito)) {
      fprintf(stderr, "bad from to values (%d %d)\n", from, to);
      return;
    }
	
    ilen = ito - ifrom + 1;

    if (bufSeq->length <= ilen) {
       bufSeq->seq = REALLOC(char, bufSeq->seq, ilen + 1);
       if (! bufSeq->seq) {
         (void) fprintf(stderr, "not enough memory\n");
         exit(10);
       }
    }
    
    (void) strncpy(bufSeq->seq, seq->seq + ifrom - 1, ilen);
    bufSeq->seq[ilen] = '\000';

    bufSeq->length = ilen;

    (void) strcpy(bufSeq->name, name ? name : seq->name);
    (void) sprintf(bufSeq->comment, "fragment %d %d %s", ifrom, ito, (rev ? "-" : "+"));
    
    if (rev) {
	(void) bio_seq_complement(bufSeq->seq);
	(void) str_reverse_string(bufSeq->seq);
    }
	
    WriteFastaSequence(stdout, bufSeq, FASTA_CHAR_PER_LINE);
}

/* ----------------------------------------------- */
/*ARGSUSED*/

main(argn, argv)
	int  argn;
	char *argv[];
{
	int nbseq, from, to, opt, ropt, iread;	
	FastaSequence *seq;
	FILE *filin;
	char strand[2], iopt[FILENAME_MAX], buffer[BUFSIZ], name[BUFSIZ];

	extern char *optarg;	/* externs for getopts (3C)	*/

	ropt = 0;
	*iopt = '\000';
	
	from = 1;
	to   = 0;
	
        while ((opt = getopt(argn, argv, "f:ht:ri:")) != -1) {

	    switch (opt) {

		case 'h':
		    (void) printf("cut fasta sequence[s]\n");
		    (void) printf("usage: cut [-f from] [-t to] [-r] [-i file]\n");
		    (void) printf(" cut fragment [from, to] (inclusive)\n");
		    (void) printf("   -f from\n");
		    (void) printf("      from > 0  : range is [from, to]\n");
		    (void) printf("      from <= 0 : range is [N-|from|, to]\n");
		    (void) printf("                N = sequence length\n");
		    (void) printf("      default: 1\n");
		    (void) printf("   -t to\n");
		    (void) printf("      to > 0  : range is [from, to]\n");
		    (void) printf("      to <= 0 : range is [from, N-|to|]\n");
		    (void) printf("                N = sequence length\n");
		    (void) printf("      default: 0\n");
		    (void) printf("   -r\n");
		    (void) printf("      reverse complement result sequence\n");
		    (void) printf("   -i file\n");
		    (void) printf("      use from,to boundaries from list file\n");
		    (void) printf("      format : from to [D|R|C|+|-  [name]] per line\n");
		    exit(0);
		    break;

		case 'f':
		    if (sscanf(optarg, "%d", &from) != 1) {
			     (void) fprintf(stderr, "bad from value\n");
			     exit(3);
			 }
			 break;

		case 't':
		    if (sscanf(optarg, "%d", &to) != 1) {
			     (void) fprintf(stderr, "bad to value\n");
			     exit(3);
			 }
			 break;
			 
		case 'r':
		    ropt = 1;
		    break;

		case 'i':
		    if (sscanf(optarg, "%s", iopt) != 1) {
			     (void) fprintf(stderr, "bad file value\n");
			     exit(3);
			 }
			 break;

		case '?':
		    (void) fprintf(stderr, "usage: cut [-f from] [-t to] [-h]\n");
		    exit(6);
		    break;
	    }
	}

	if (*iopt && (! (filin = fopen(iopt, "r")))) {
	   (void) fprintf(stderr, "%s: file not found\n", iopt);
	   exit(7);
	}
	
	seq = NewFastaSequence();
	
	nbseq = 0;
	
	if (! *iopt) {
	
	
		while (ReadFastaSequence(stdin, seq)) {
	
		    nbseq++;
	
		    if (! seq->ok) {
			(void) fprintf(stderr, "error at seq # %d\n", nbseq);
			continue;
		    }
		    
		    sPrintSequence(seq, from, to, ropt, seq->name);
		}
	
	}
	
	else {
		    ReadFastaSequence(stdin, seq);
		    
		    nbseq++;
	
		    if (! seq->ok) {
			(void) fprintf(stderr, "error at seq # %d\n", nbseq);
			exit(8);
		    }
		    
		    while (fgets(buffer, sizeof(buffer), filin)) {
		       
		      iread = sscanf(buffer, "%d%d%1s%s", &from, &to, strand, name);
		      if (iread < 2) {
		        fprintf(stderr, "ignored boundaries at line \"%s\"\n", buffer);
		        continue;
		      }
		      if (iread >= 3) {
		        ropt = (((*strand == '-') || (*strand == 'R') || (*strand == 'C')) ? 1 : 0);
		      }
		      if (iread < 4) {
		        (void) strcpy(name, seq->name);
		      }
		      
		      sPrintSequence(seq, from, to, ropt, name);
		    }
		    
		    fclose(filin);
	}

	FreeFastaSequence(seq);
	
	exit(0);
	
	/*NOTREACHED*/
}


	


