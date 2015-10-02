
#ifndef FASTA_SEQ_WRITER_H_
#define FASTA_SEQ_WRITER_H_

#include "sequence.h"


void printOnlySeqFromFastaSeqPtr(fastaSeqPtr, FILE*);

void printOnlySeqFromChar(char*, FILE*);

void printOnlyHeaderFromFastaSeqPtr(fastaSeqPtr, FILE*);

void printOnlyHeaderFromTable(element_from_header*, FILE*);

void printHeaderAndSeqFromFastaSeqPtr(fastaSeqPtr, FILE*);


#endif
