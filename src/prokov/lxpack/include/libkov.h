/* ---------------------------------------------------------------- */
/* Copyright (c) Atelier de BioInformatique                         */
/* @file: libkov.h                                                  */
/* @desc: libkov - include file                                     */
/*                                                                  */
/* @history:                                                        */
/* @+   <Fred Nikitin + Marc Heuveline> : Aug 99 : first version    */
/* @+   <Gloup> : Oct 99 : last revised version                     */
/* ---------------------------------------------------------------- */

#ifndef _H_libkov

#define _H_libkov

#define K_MAX      8            /* kuple maxi                   */
#define POW_KMAX   65536        /*  = 4^KMAX                    */
#define POW_KMIN   16384        /*  = 4^(KMAX-1)                */

#define MAX_WIN    4096         /* taille maxi fenetre          */

                                /* --- default values ----      */
#define KUP_DFT    3            /* default kuple                */
#define PCODANT    80           /* proba codant A Priori (%)    */
#define WIN_WIDTH  120          /* window width                 */
#define WIN_STEP   30           /* window step                  */
#define MIN_CDS    120          /* minimun cds length           */
#define THRESH_CDS 40           /* proba cds min (%)            */
#define MATNAME    "Matrix"     /* nom fichier de sortie        */

#define START_DFT  "/ATG/GTG/TTG/"
#define STOP_DFT   "/TAA/TAG/TGA/"

#define VERSION    "v1.2"       /* version number               */
#define VERSIZE    16           /* header version size          */

                                /* env. variable for Matrix dir */
#define PROKOV_MAT_ENV "PROKOV_MATDIR"

/* ==================================================== */
/* structures                                           */
/* ==================================================== */

typedef struct {
        int     kupleC, kupleN;
        int     powkC, powkN;

        int     occpos[3][POW_KMAX],    /* occurences   */
                occneg[3][POW_KMAX],
                occnc[POW_KMAX];

        float   probpos[3][POW_KMAX],   /* probas cond. */
                probneg[3][POW_KMAX],
                probnc[POW_KMAX];
                
        float   pinipos[3][POW_KMIN],   /* probas init. */
                pinineg[3][POW_KMIN],
                pininc[POW_KMIN];

        char    version[VERSIZE];
} MarkovMatrix;

typedef struct {
        float   probpos[3],
                probneg[3], 
                probnc;
} ProbaArray;

        
/* ==================================================== */
/* macros                                               */
/* ==================================================== */

#define MIN(x, y) (((x) < (y)) ? (x) :  (y))
#define MAX(x, y) (((x) > (y)) ? (x) :  (y))
#define ABS(x)    (((x) > 0)   ? (x) : -(x))

#ifndef NEW
#define NEW(typ)                (typ*)malloc(sizeof(typ)) 
#define NEWN(typ, dim)          (typ*)malloc((unsigned)(dim) * sizeof(typ))
#define REALLOC(typ, ptr, dim)  (typ*)realloc((void *) (ptr), (unsigned)(dim) * sizeof(typ))
#define FREE(ptr)               free(ptr)
#endif

/* ==================================================== */
/* prototypes                                           */
/* ==================================================== */

                        /* libkov_util.c                */

int     Erreur           (char *msg, int stat);
char    *Complement      (char *str);
char    *Reverse         (char *str);
char    *Upper           (char *str);
char    *RNAtoDNA        (char *str);
int     Valeur           (char s);
int     ValeurComplement (int valeur);
int     Pow4             (int n);
int     Codage           (char *seq, int kuple);
int     CodageComplement (int index, int kuple);
int     SuiteCodage      (char *seq, int prec, int kuple);
int     FindCodon        (char *seq, int lenseq, char *codons, int from, int to);
int     InternalFrame   (int dframe, int direct, int seqlen);

void    ComputeNegOccurences(int occneg[3][POW_KMAX],
                             int occpos[3][POW_KMAX], int kuple);


                        /* libkov_io.c                  */

int     AssignToStdin   (char *filename);
FILE    *OpenFile       (char *filename, char *mode);
char    *MatrixPathName (char *shortname);

int     WriteMatrix     (MarkovMatrix *mat, int isAscii, FILE *stream);
int     ReadMatrix      (MarkovMatrix *mat, int isAscii, FILE *stream);

                        /* libkov_proba.c               */

void    ProbaCond   (MarkovMatrix *mat);
                  
void    ProbaMarkov (char *seq, int seqlen, MarkovMatrix *mat,
                     ProbaArray *win);
                     
void    ProbaBayes  (ProbaArray *win, float prior, ProbaArray *bay);

#endif
