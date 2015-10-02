/* ==================================================== */
/* @File kimono.h                                       */
/* @history:                                            */
/* @+   Nov. 03 <Gloup> first draft                     */
/* ==================================================== */

#define _H_kimono

#ifndef _H_fasta_io
#include "fasta_io.h"
#endif

/* -------------------------------------------- */
/* constantes                                   */
/* -------------------------------------------- */

#define NB_NUC          4
#define NB_AA	        21
#define NB_CODONS       64

#define SHIFT_NUC       2        /* # bits for NB_NUC */

#define KIM_DNA 	"ACGT"
#define KIM_AA 		"ACDEFGHIKLMNPQRSTVWY*"

#define MASKED_SYMBOL   'X'

/* -------------------------------------------- */
/* struct de donnees                            */
/* -------------------------------------------- */

typedef struct {
  int from, to, strand;
} Interval;

typedef struct {
  int size, capacity;
  Interval *interval;
} SetOfIntervals;

typedef struct {
  char index; 
  char mask;
} SeqCell;

typedef struct {
  int length;
  char symbol[128];
} Alphabet;

typedef struct {
  Alphabet *alphabet;
  FastaSequence *fasta;
  SeqCell *cell;
} Sequence;

/* -------------------------------------------- */
/* counting, indexing and distance functions    */
/* -------------------------------------------- */

typedef void  (*IndexFunction) (Sequence *seq, SetOfIntervals *zones, int phase);
typedef void  (*CountFunction) (Sequence *seq, int from, int to, int count[][2], int nbins, int increment);
typedef void  (*MixFunction)   (int left[][2], int right[][2], int countLeft[], int countRight[], int nbins);
typedef float (*DistFunction)  (int bins1[], int bins2[], int nbins, int ibin);

/* -------------------------------------------- */
/* masks for sequence                           */
/* -------------------------------------------- */

#define EXCLUDE_MASK   0x01     /* 00000001				 */
#define INCLUDE_MASK   0x7e     /* 01111110				 */
#define DIRECT_MASK    0x02     /* 00000010				 */
#define REVERSE_MASK   0x04     /* 00000100				 */

/* -------------------------------------------- */
/* macros                                       */
/* -------------------------------------------- */

#define IS_UPPER(c)         (((c) >= 'A') && ((c) <= 'Z'))
#define IS_LOWER(c)         (((c) >= 'a') && ((c) <= 'z'))
#define IS_ALPHA(c)         (IS_UPPER(c) || IS_LOWER(c))

#define TO_UPPER(c)         ((c) - 'a' + 'A')
#define TO_LOWER(c)         ((c) - 'A' + 'a')

#define UPPER(x)   	    (IS_LOWER(x) ? TO_UPPER(x) : (x))
#define LOWER(x)   	    (IS_UPPER(x) ? TO_LOWER(x) : (x))

#define MIN2(a, b)          ((a) < (b) ? (a) : (b))     
#define MAX2(a, b)          ((a) > (b) ? (a) : (b))     

#ifndef NEW
#define NEW(type)               (type *) malloc(sizeof(type))
#define NEWN(type, n)           (type *) malloc((n) * sizeof(type))
#define REALLOC(type, ptr, n)   (type *) realloc(ptr, (n) * sizeof(type))
#define FREE(ptr)               free(ptr)
#endif

/* -------------------------------------------- */
/* Prototypes                                   */
/* -------------------------------------------- */

                                  /* kim_util.c */
                                
int     Erreur          (char *msg , int stat);
int     MemoryErreur    (char *where, int stat);
int     AssignToStdin   (char *filename);

                                  /* kim_alphabet.c */
                                  
Alphabet *FreeAlphabet    (Alphabet *alpha);
Alphabet *NewAlphabet     (char *symbols);
Alphabet *CopyAlphabet    (Alphabet *alpha);
Alphabet *NucleicAlphabet ();
Alphabet *AminoAlphabet   ();

                                  /* kim_sequence.c */

Sequence *FreeSequence   (Sequence *seq);
Sequence *NewSequence    (Alphabet *alpha);
int      ReadSequence    (FILE *streamin, Sequence *seq);

                                 /* kim_genetic.c */
                                 
int 	BaseComplement	(int ch);
char 	*SeqComplement	(char *str);
int 	CodonTranslate	(char *codon, int codid);

                                  /* kim_zone.c */
                                  
SetOfIntervals *FreeZones  (SetOfIntervals *set);
SetOfIntervals *ReadZones  (int seqlength);

                                  /* kim_indexer.c */

void NoZoneIndexer  (Sequence *seq, SetOfIntervals *zones, int phase);
void ExcludeIndexer (Sequence *seq, SetOfIntervals *zones, int phase);
void GeneIndexer    (Sequence *seq, SetOfIntervals *zones, int phase);
void PhaseIndexer   (Sequence *seq, SetOfIntervals *zones, int phase);
void CodonIndexer   (Sequence *seq, SetOfIntervals *zones, int phase);

                                  /* kim_counter.c */

void BasicCounter      (Sequence *seq, int from, int to, 
                        int count[][2], int nbins, int increment);
void CodonSkewCounter  (Sequence *seq, int from, int to,
                        int count[][2], int nbins, int increment);

                                  /* kim_mixer.c */
void BasicMixer        (int left[][2], int right[][2], 
                        int *countLeft, int *countRight, int nbins);
void CrossMixer        (int left[][2], int right[][2], 
                        int *countLeft, int *countRight, int nbins);

                                  /* kim_distance.c */
                                  
float 	ChisTwo(int bins1[], int bins2[], int nbins, int ibin);
float 	Euclid(int bins1[], int bins2[], int nbins, int ibin);
float 	HellTwo(int bins1[], int bins2[], int nbins, int ibin);
int     SignTwoBins(int bins1[], int bins2[]);

				/* kim_codonskew.c */
				
void	CodonSkew	(int code, float score[NB_NUC][NB_CODONS]);

                                /* kim*_help.c */
void    PrintHelp (void);
int     ExitUsage (int stat);

