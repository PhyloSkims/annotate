#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "sequence.h"
#include "fasta_header_parser.h"


void printOnlySeqFromFastaSeqPtr(fastaSeqPtr seq, FILE* output)
{
	char nuc;
	int n=60;
	int l = strlen(seq->sequence);
	for (n=60; n<l; n+=60)
	{
		nuc = seq->sequence[n];
		seq->sequence[n]=0;
		fprintf(output,"%s\n",seq->sequence+n-60);
		seq->sequence[n]=nuc;
	}
	fprintf(output,"%s\n",seq->sequence+n-60);
}


void printOnlySeqFromChar(char* seq, FILE* output)
{
	char nuc;
	int n=60;
	int l = strlen(seq);
	for (n=60; n<l; n+=60)
	{
		nuc = seq[n];
		seq[n]=0;
		fprintf(output,"%s\n",seq+n-60);
		seq[n]=nuc;
	}
	fprintf(output,"%s\n",seq+n-60);
}


void printOnlyHeaderFromFastaSeqPtr(fastaSeqPtr seq, FILE* output)
{
	fprintf(output,">%s\n",seq->rawheader);
}


void printOnlyHeaderFromTable(element_from_header* header, FILE* output)
{
	int i;
	int nbf;

	nbf = atoi(header[0].value);

	fprintf(output,">%s ",header[1].value);

	for (i = 2; i <= nbf; i++)
	{
		if (strcmp(header[i].name, "definition") != 0)
		{
			fprintf(output,"%s",header[i].name);
			fprintf(output,"=");
		    fprintf(output,"%s; ",header[i].value);
		}
	}

	if (strcmp(header[nbf].name, "definition") == 0)
		fprintf(output,"%s; ",header[nbf].value);

	fprintf(output,"\n");
}


void printHeaderAndSeqFromFastaSeqPtr(fastaSeqPtr seq, FILE* output)
{
	printOnlyHeaderFromFastaSeqPtr(seq, output);
	printOnlySeqFromFastaSeqPtr(seq, output);
}
