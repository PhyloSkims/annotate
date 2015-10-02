/* ---------------- */
/* @file: ctest.c   */
/* ---------------- */

#include <stdio.h>
#include <stdlib.h>

#include "libbio.h"

#define BUFFER_SIZE	 1024
#define CPL		 10

#define PRINT(s)	 sPrintSequence(stdout, (s), CPL);

/* ---------------------------- */

static void sPrintSequence(FILE *fd, char *seq, int nchars)
{
	long count;

	for (count = 0 ; *seq ; seq++, count++)  {
	    if (count && ((count % nchars) == 0))
		fputc('\n', fd);
	    fputc(*seq, fd);
	}

	fputc('\n', fd);
}

/* ---------------------------- */

main ()
{
	long len;
	char line[256], seq[BUFFER_SIZE], sseq[BUFFER_SIZE];


	*seq = '\000';
		
	while (gets(line)) 
	    strcat(seq, line);

	len = strlen(seq);

	printf ("[read] %ld symbols\n", len);
	PRINT(seq);

	printf("[pwg_trim_leading]\n");
	strcpy(sseq, "    ");
	strcat(sseq, seq);
	PRINT(sseq);
	str_trim_leading(sseq);
	PRINT(sseq);

	printf("[str_trim_trailing]\n");
	strcpy(sseq, seq);
	strcat(sseq,  "   ");
	PRINT(sseq);
	str_trim_trailing(sseq);
	PRINT(sseq);

	printf("[str_pad_right]\n");
	strcpy(sseq, seq);
	str_pad_right(sseq, len + len/2, 'x');
	PRINT(sseq);

	printf("[str_pad_left]\n");
	strcpy(sseq, seq);
	str_pad_left(sseq, len + len/2, 'x');
	PRINT(sseq);

	printf("[str_drop_string]\n");
	strcpy(sseq, seq);
	str_drop_string(sseq, 10, 5);
	PRINT(sseq);

	printf("[str_insert_string]\n");
	strcpy(sseq, seq);
	str_insert_string(sseq, "xxxxx", 10);
	PRINT(sseq);

	printf("[str_extract_string]\n");
	str_extract_string(sseq, seq, 10, 20);
	PRINT(sseq);


	printf("[str_extract_to_mark]\n");
	str_extract_to_mark(sseq, seq, 10, 'C');
	PRINT(sseq);

	printf("[bio_seq_complement]\n");
	strcpy(sseq, seq);
	bio_seq_complement (sseq);
	PRINT(sseq);

	printf("[str_reverse_string]\n");
	strcpy(sseq, seq);
	str_reverse_string(sseq);
	PRINT(sseq);

	printf("[end of test]\n");

	exit(0);

}
