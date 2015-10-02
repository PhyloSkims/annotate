/* ---------------------------------------------------------------- */
/* @file: libkov_io.c                                               */
/* @desc: Prokov Input/Output                                       */
/*                                                                  */
/* @history:                                                        */
/* @+   <ABI> : Aug 99 : first version                              */
/* @+   <Gloup> : Oct 99 : last revised version                     */
/* @+   <Gloup> : Dec 00 : HelixWare port - V1.2                    */
/* ---------------------------------------------------------------- */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "libkov.h"

#ifndef FILENAME_MAX
#define FILENAME_MAX 1024
#endif

#ifdef LITTLE_ENDIAN
#define ENDIAN(n) sEndian(n)
#else
#define ENDIAN(n) n
#endif

#define LONG_OUTPUT 0   /* should be 0 for compatibility with V1.2 and earlier */

/* ============================================ */
/* Machine binary conversions                   */
/* ============================================ */

/* -------------------------------------------- */
/* static                                       */
/* convert int  little <-> big endian           */
/* binary files are coded in big-endian format  */
/* -------------------------------------------- */
int sEndian (int n)
{
    union {
       int  integer;
       unsigned char byte[4];
    } convert;
    
    unsigned char tempo;

    convert.integer = n;
    
    tempo           = convert.byte[0];
    convert.byte[0] = convert.byte[3];
    convert.byte[3] = tempo;
    tempo           = convert.byte[1];
    convert.byte[1] = convert.byte[2];
    convert.byte[2] = tempo;
    
    return convert.integer;
}

/* -------------------------------------------- */
/* static                                       */
/* convert array src to dst big endian          */
/* -------------------------------------------- */
static void sArrayEndian(int dst[], int src[], int taille)
{
    int i;
    
    for (i = 0 ; i < taille ; i++)
        dst[i] = ENDIAN(src[i]);
}

/* ============================================ */
/* File Utilities                               */
/* ============================================ */

/* -------------------------------------------- */
/* open a file and assign it to stdin           */
/* -------------------------------------------- */

int AssignToStdin(char *filename)
{
    char buffer[256];
        
    if (! freopen(filename, "r", stdin)) {
        (void) sprintf(buffer, "cant open file %s to read", filename);
        return Erreur(buffer, 0);
    }

    return 1;
}

/* -------------------------------------------- */
/* open a file                                  */
/* -------------------------------------------- */

FILE *OpenFile(char *filename, char *mode)
{
    FILE *stream;
    char buffer[256];
        
    if (! (stream = fopen(filename, mode))) {
        (void) sprintf(buffer, "cant open file %s in mode %s",
                        filename, mode);
        (void) Erreur(buffer, 0);
    }

    return stream;
}

/* -------------------------------------------- */
/* get fullpathname for matrix file             */
/* -------------------------------------------- */

char * MatrixPathName(char *shortname) {

    char *env;
    
    static char sPath[FILENAME_MAX];

    env = getenv(PROKOV_MAT_ENV);
    
    if (env)
        (void) strcat(strcpy(sPath, env), "/");
    else
        (void) strcpy(sPath, "");

    return strcat(sPath, shortname);
}


/* ============================================ */
/* Matrix io                                    */
/* ============================================ */

/* -------------------------------------------- */
/* static                                       */
/* ensure count of occurences is not 0          */
/* -------------------------------------------- */
static void sCheckOccurences(int occ[], int taille)
{
    int i;
    
    for (i = 0 ; i < taille ; i++) {
        if (occ[i] == 0)
            occ[i] = 1;
    }
}

/* -------------------------------------------- */
/* write Header - Ascii format                  */
/* -------------------------------------------- */

static int sWriteMatrixHeaderAscii(MarkovMatrix *mat, FILE *stream)
{
    int ok;
    
    ok = fprintf(stream, "%s %d %d %d %d\n",
                 mat->version, mat->kupleC, mat->kupleN,
                 mat->powkC, mat->powkN);

    return (ok > 0);
}

/* -------------------------------------------- */
/* write Header - Binary format                 */
/* -------------------------------------------- */

static int sWriteMatrixHeaderBinary(MarkovMatrix *mat, FILE *stream)
{
    int  ok, nb, val;
    char buffer[VERSIZE];

    (void) memset(buffer, 0, VERSIZE);
    (void) strcpy(buffer, mat->version);  /* version : 16 bytes */    
    nb = fwrite(buffer, sizeof(char), VERSIZE, stream);
    ok = (nb == VERSIZE); 

    val = ENDIAN(mat->kupleC);            /* kupleC : 4 bytes   */
    nb = fwrite(&val, sizeof(int), 1, stream);
    ok = ok && (nb == 1); 

    val = ENDIAN(mat->kupleN);            /* kupleN : 4 bytes   */
    nb = fwrite(&val, sizeof(int), 1, stream);
    ok = ok && (nb == 1); 

    val = ENDIAN(mat->powkC);             /* powkC : 4 bytes    */
    nb = fwrite(&val, sizeof(int), 1, stream);
    ok = ok && (nb == 1); 

    val = ENDIAN(mat->powkN);            /* powkN : 4 bytes     */
    nb = fwrite(&val, sizeof(int), 1, stream);
    ok = ok && (nb == 1); 

    return ok;
}


/* -------------------------------------------- */
/* write Matrix Body - Ascii format             */
/* -------------------------------------------- */

static int sWriteMatrixBodyAscii(MarkovMatrix *mat, FILE *stream)
{
    int i, phase, ok;
    
    ok = 1;
    
    for (i = 0 ; i < mat->powkC ; i++) {
       
       for (phase = 0 ; phase < 3 ; phase++) {
           ok = ok && (fprintf(stream, "%d\t", mat->occpos[phase][i]) > 0);

#if LONG_OUTPUT
           ok = ok && (fprintf(stream, "%d\t", mat->occneg[phase][i]) > 0);
#endif
       }

       (void) fprintf(stream, "\n");
    }

    for (i = 0 ; i < mat->powkN ; i++)
        ok = ok && (fprintf(stream, "%d\n", mat->occnc[i]) > 0);

   return ok;
}

/* -------------------------------------------- */
/* write Matrix Body - Binary format            */
/* -------------------------------------------- */

static int sWriteMatrixBodyBinary(MarkovMatrix *mat, FILE *stream)
{
    int phase, nb;
    int buffer[POW_KMAX];
    
    nb    = 0;
    
    for (phase = 0 ; phase < 3 ; phase++) {
        sArrayEndian(buffer, mat->occpos[phase], mat->powkC);
        nb += fwrite(buffer, sizeof(int), mat->powkC, stream);

#if LONG_OUTPUT
        sArrayEndian(buffer, mat->occneg[phase], mat->powkC);
        nb += fwrite(buffer, sizeof(int), mat->powkC, stream);
#endif
    }

    sArrayEndian(buffer, mat->occnc, mat->powkN);
    nb += fwrite(buffer, sizeof(int), mat->powkN, stream);

#if LONG_OUTPUT
    return (nb == ((mat->powkC * 6) + mat->powkN));
#else
    return (nb == ((mat->powkC * 3) + mat->powkN));
#endif
}

/* -------------------------------------------- */
/* read Header - Ascii format                   */
/* -------------------------------------------- */

static int sReadMatrixHeaderAscii(MarkovMatrix *mat, FILE *stream)
{
    int  ok;
    
    ok = fscanf(stream, "%s%d%d%d%d\n",
                mat->version, &(mat->kupleC), &(mat->kupleN), 
                              &(mat->powkC),  &(mat->powkN));

    ok =    (ok > 0) 
         && (mat->kupleC <= K_MAX) && (mat->kupleN <= K_MAX)
         && (mat->kupleC > 0)      && (mat->kupleN > 0)
         && (mat->powkC == Pow4(mat->kupleC)) 
         && (mat->powkN == Pow4(mat->kupleN))
         && (! strcmp(mat->version, VERSION));
         
    return ok;
}

/* -------------------------------------------- */
/* read Header - Binary format                  */
/* -------------------------------------------- */

static int sReadMatrixHeaderBinary(MarkovMatrix *mat, FILE *stream)
{
    int  ok, nb, val;

                                        /* version : 16 bytes   */    
    nb = fread(mat->version, sizeof(char), VERSIZE, stream);
    ok =    (nb == VERSIZE)
         && (! strcmp(mat->version, VERSION));

                                        /* kupleC : 4 bytes     */
    nb = fread(&val, sizeof(int), 1, stream);
    mat->kupleC = ENDIAN(val);
    ok =    ok 
         && (nb == 1)
         && (mat->kupleC <= K_MAX)
         && (mat->kupleC > 0);

                                        /* kupleN : 4 bytes     */
    nb = fread(&val, sizeof(int), 1, stream);
    mat->kupleN = ENDIAN(val);
    ok =    ok 
         && (nb == 1)
         && (mat->kupleN <= K_MAX)
         && (mat->kupleN > 0);

                                        /* powkC : 4 bytes      */
    nb = fread(&val, sizeof(int), 1, stream);
    mat->powkC = ENDIAN(val);
    ok =    ok 
         && (mat->powkC = Pow4(mat->kupleC));

                                        /* powkC : 4 bytes      */
    nb = fread(&val, sizeof(int), 1, stream);
    mat->powkN = ENDIAN(val);
    ok =    ok 
         && (mat->powkN = Pow4(mat->kupleN));

    return ok;
}

/* -------------------------------------------- */
/* read Matrix Body - Ascii format              */
/* -------------------------------------------- */

static int sReadMatrixBodyAscii(MarkovMatrix *mat, FILE *stream)
{
    int i, phase, ok;
    
    ok = 1;
    
    for (i = 0 ; i < mat->powkC ; i++) {
       
       for (phase = 0 ; phase < 3 ; phase++) {
           ok = ok && (fscanf(stream, "%d\t", &(mat->occpos[phase][i])) > 0);

#if LONG_OUTPUT
           ok = ok && (fscanf(stream, "%d\t", &(mat->occneg[phase][i])) > 0);
#endif
       }
    }

#if ! LONG_OUTPUT
    ComputeNegOccurences(mat->occneg, mat->occpos, mat->kupleC);
#endif
   
    for (i = 0 ; i < mat->powkN ; i++)
       ok = ok && (fscanf(stream, "%d\n", &(mat->occnc[i])) > 0);
   
   if (ok) {
        for (phase = 0 ; phase < 3 ; phase++)
            sCheckOccurences(mat->occpos[phase], mat->powkC);
        for (phase = 0 ; phase < 3 ; phase++)
            sCheckOccurences(mat->occneg[phase], mat->powkC);
        sCheckOccurences(mat->occnc, mat->powkN);
   }
   
   return ok;
}

/* -------------------------------------------- */
/* read Matrix Body - Binary format             */
/* -------------------------------------------- */

static int sReadMatrixBodyBinary(MarkovMatrix *mat, FILE *stream)
{
    int phase, nb, ok;
    int buffer[POW_KMAX];
    
    nb = 0;
    
    for (phase = 0 ; phase < 3 ; phase++) {
        nb += fread(buffer, sizeof(int), mat->powkC, stream);
        sArrayEndian(mat->occpos[phase], buffer, mat->powkC);

#if LONG_OUTPUT
        nb += fread(buffer, sizeof(int), mat->powkC, stream);
        sArrayEndian(mat->occneg[phase], buffer, mat->powkC);
#endif
    }

#if ! LONG_OUTPUT
    ComputeNegOccurences(mat->occneg, mat->occpos, mat->kupleC);
#endif

    nb += fread(buffer, sizeof(int), mat->powkN, stream);
    sArrayEndian(mat->occnc, buffer, mat->powkN);

#if LONG_OUTPUT
    ok = (nb == ((mat->powkC * 6) + mat->powkN));
#else
    ok = (nb == ((mat->powkC * 3) + mat->powkN));
#endif
    
    if (ok) {
        for (phase = 0 ; phase < 3 ; phase++)
            sCheckOccurences(mat->occpos[phase], mat->powkC);
        for (phase = 0 ; phase < 3 ; phase++)
            sCheckOccurences(mat->occneg[phase], mat->powkC);
        sCheckOccurences(mat->occnc, mat->powkN);
    }

    return ok;
}

/* -------------------------------------------- */
/* write/read Matrix - Interface calls          */
/* -------------------------------------------- */

int WriteMatrix(MarkovMatrix *mat, int isAscii, FILE *stream)
{

    if (isAscii) {
       if (! sWriteMatrixHeaderAscii(mat, stream)) {
          Erreur("Invalid ascii matrix header", 0);
          return 0;
       }
       if (! sWriteMatrixBodyAscii(mat, stream)) {
          Erreur("Invalid ascii matrix body", 0);
          return 0;
       }
    }
    else {
       if (! sWriteMatrixHeaderBinary(mat, stream)) {
          Erreur("Invalid binary matrix header", 0);
          return 0;
       }
       if (! sWriteMatrixBodyBinary(mat, stream)) {
          Erreur("Invalid binary matrix body", 0);
          return 0;
       }
    }

    return 1;
}

int ReadMatrix(MarkovMatrix *mat, int isAscii, FILE *stream)
{

    if (isAscii) {
       if (! sReadMatrixHeaderAscii(mat, stream)) {
          Erreur("Invalid ascii matrix header", 0);
          return 0;
       }
       if (! sReadMatrixBodyAscii(mat, stream)) {
          Erreur("Invalid ascii matrix body", 0);
          return 0;
       }
    }
    else {
       if (! sReadMatrixHeaderBinary(mat, stream)) {
          Erreur("Invalid binary matrix header", 0);
          return 0;
       }
       if (! sReadMatrixBodyBinary(mat, stream)) {
          Erreur("Invalid binary matrix body", 0);
          return 0;
       }
    }

    return 1;
}

