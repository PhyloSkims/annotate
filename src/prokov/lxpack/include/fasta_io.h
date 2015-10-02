/* ---------------------------------------------------------------- */
/* @file: fasta_io.h                                                */
/* @desc: sequence IO in fasta format / include file                */
/*                                                                  */
/* @history:                                                        */
/* @+       <Gloup> : Aug 92 : first version                        */
/* @+       <Gloup> : Nov 95 : last revised version                 */
/* ---------------------------------------------------------------- */

#ifndef _H_fasta_io

#define _H_fasta_io

/* ==================================================== */
/* Constantes                                           */
/* ==================================================== */

#define FASTA_NAMLEN  64        /* max length of seq. name       */
#define FASTA_COMLEN  512       /* max length of seq. comment    */

#define FASTA_CHAR_PER_LINE 50  /* # of chars per line in output */

/* ==================================================== */
/* Macros standards                                     */
/* ==================================================== */

#ifndef NEW
#define NEW(typ)                (typ*)malloc(sizeof(typ)) 
#define NEWN(typ, dim)          (typ*)malloc((unsigned)(dim) * sizeof(typ))
#define REALLOC(typ, ptr, dim)  (typ*)realloc((void *) (ptr), (unsigned long)(dim) * sizeof(typ))
#define FREE(ptr)               free(ptr)
#endif

/* ==================================================== */
/* Structures de donnees                                */
/* ==================================================== */

typedef struct {                        /* -- Sequence ---------------- */
        int     ok;                     /* error flag                   */
        long    length,                 /* longueur                     */
                offset,                 /* offset                       */
                bufsize;                /* size of current seq buffer   */
        char    name[FASTA_NAMLEN],     /* nom                          */
                comment[FASTA_COMLEN],  /* commentaire                  */
                *seq;                   /* sequence                     */
} FastaSequence, *FastaSequencePtr;

/* ==================================================== */
/*  Prototypes                                          */
/* ==================================================== */

                                        /* libfasta.c   */

char             *GetFastaName      ( char *buffer );
char             *GetFastaComment   ( char *buffer );

FastaSequencePtr NewFastaSequence   ( void );
FastaSequencePtr FreeFastaSequence  ( FastaSequencePtr seq );

int              ReadFastaSequence  ( FILE *streamin, FastaSequencePtr seq );
int              GetFastaSequence   ( FILE *streamin, FastaSequencePtr seq );
void             WriteFastaSequence ( FILE *streamou, FastaSequencePtr seq , int char_per_line );

void             RewindFastaDB      ( FILE *streamin );

#endif
