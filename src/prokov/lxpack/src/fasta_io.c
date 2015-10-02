/* ---------------------------------------------------------------- */
/* @file: fasta_io.c                                                */
/* @desc: sequence IO in fasta format                               */
/*                                                                  */
/* @history:                                                        */
/* @+       <Gloup> : Aug 92 : first version                        */
/* @+       <Gloup> : Nov 95 : last revised version                 */
/* ---------------------------------------------------------------- */

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "fasta_io.h"

#define SECURE   0      /* 1 : secure copy ; 0 : quickest copy  */

#define READ_NEXT 0
#define PUSH_BACK 1

#define SERIAL    1
#define INDEXED   0

#ifdef MAC_OS_C
#define LINE_FEED '\r'
#else
#define LINE_FEED '\n'
#endif

static int sRetained = 0;

/* -------------------------------------------- */
/* @static: lecture bufferisee                  */
/* -------------------------------------------- */
static char * sNextIOBuffer(FILE *streamin, int retain, int serial)
{
        char  *buf, *end;
        int   reused;
        
        static char sBuffer[BUFSIZ];    /* BUFSIZ in <stdio.h>  */
        
        reused = (retain || sRetained) && serial;

        sRetained = retain;

        buf = (   reused 
                ? sBuffer 
                : fgets(sBuffer, sizeof(sBuffer), streamin));           

        if (! buf)                /* end of file            */
            return NULL;

        end = buf + strlen(buf) - 1;

        if (*end == LINE_FEED)    /* remove trailing \n     */
           *end = '\000';

        return (  *buf || reused
                ? buf
                : sNextIOBuffer(streamin, retain, serial));
}

#if SECURE
/* -------------------------------------------- */
/* compte le nombre de caracteres alpha dans    */
/* un buffer                                    */
/* -------------------------------------------- */
static long sCountAlpha(char *buf)
{
        long count;
        
        for (count = 0 ; *buf ; buf++)
            if (isalpha(*buf))
                count++;
        
        return count;
}
#endif

#if SECURE
/* -------------------------------------------- */
/* copy only alpha chars from s2 to s1          */
/* -------------------------------------------- */
static char *sStrcpyAlpha(char *s1, char *s2)
{
        for( ; *s2 ; s2++)
            if (isalpha(*s2))
                *s1++ = *s2;

        *s1 = '\000';

        return s1;
}
#endif

/* -------------------------------------------- */
/* skip to next space in buffer                 */
/* -------------------------------------------- */
static char *sNextSpace(char *buffer)
{
        for (; *buffer ; buffer++)
           if (isspace(*buffer))
                return buffer;
        
        return NULL;
}

/* -------------------------------------------- */
/* returns sequence name (FASTA)                */
/* -------------------------------------------- */
char *GetFastaName(char *buffer)
{
        char c[2];
        static char name[FASTA_NAMLEN];

        *c = buffer[FASTA_NAMLEN];
         buffer[FASTA_NAMLEN] = '\000';
        
        if (sscanf(buffer + 1, "%s", name) != 1)
            (void) strcpy(name, "<no Name>");

        buffer[FASTA_NAMLEN] = *c;
         
        return name;
}

/* -------------------------------------------- */
/* returns sequence comment (FASTA)             */
/* -------------------------------------------- */
char *GetFastaComment(char *buffer)
{
        char   *space;
        static char comment[FASTA_COMLEN];

        buffer[FASTA_COMLEN] = '\000';
        
        space = sNextSpace(buffer);
        
        (void) strcpy(comment, (space ? space + 1 : "\000"));

        return comment;
}

/* -------------------------------------------- */
/* liberation d'une sequence                    */
/* -------------------------------------------- */
FastaSequencePtr FreeFastaSequence(FastaSequencePtr seq)
{
        if (seq) {
            if (seq->seq)  FREE(seq->seq);
            FREE(seq);
        }

        return NULL;
}
        
/* -------------------------------------------- */
/* allocation d'une sequence                    */
/* -------------------------------------------- */
FastaSequencePtr NewFastaSequence()
{
        FastaSequencePtr seq;
        
        if (! (seq = NEW(FastaSequence)))
            return NULL;
           
        seq->length   = 0;
        seq->offset   = -1;

        if (! (seq->seq = NEWN(char,  BUFSIZ)))
            return FreeFastaSequence(seq);

        seq->bufsize = BUFSIZ;

        *(seq->name)    = '\000';
        *(seq->comment) = '\000';

        seq->ok = 1;
        
        return seq;
}

/* -------------------------------------------- */
/* lecture/redimensionnement d'une sequence au  */
/* format Fasta                                 */
/* <Lecture Serie>                              */
/* returns : 0 -> last sequence                 */
/*           1 -> more to read                  */
/*           <but> you must check seq->ok       */
/*           that may be '0' (memory error)     */
/* -------------------------------------------- */
int ReadFastaSequence(FILE *streamin, FastaSequencePtr seq)
{
        long    readlen, buflen;
        char    *buffer, *tbuf;

        seq->ok = 0;                             /* assume error         */

        buflen = seq->length = 0; 
        
        seq->offset = ftell(streamin);

        buffer = sNextIOBuffer(streamin, READ_NEXT, SERIAL);

        if (! (buffer && (*buffer == '>')))     /* sync error           */
            return 0;                           /* last sequence        */
        
        if (seq->offset)
            seq->offset -= (strlen(buffer) + 1);

        (void) strcpy(seq->name,    GetFastaName(buffer));
        
        (void) strcpy(seq->comment, GetFastaComment(buffer));
        
        while (buffer = sNextIOBuffer(streamin, READ_NEXT, SERIAL)) {

            if (*buffer == '>') {
                (void) sNextIOBuffer(streamin, PUSH_BACK, SERIAL); /* push it back */
                break;
            }

#if SECURE          
            readlen = sCountAlpha(buffer);
#else
            readlen = strlen(buffer);
#endif
        
            buflen +=  readlen;
   
            if (buflen >= seq->bufsize) {
            
                if (! (tbuf = REALLOC(char, seq->seq, 2 * buflen + 1)))
                   return 1;  /* but seq->ok is 0  */

                seq->seq = tbuf;
                
                seq->bufsize = 2 * buflen + 1;
                
            }           
#if SECURE
            sStrcpyAlpha(seq->seq + seq->length, buffer);
#else
            (void) memcpy(seq->seq + seq->length, buffer, readlen);
#endif
        
            seq->length = buflen;
        
        }

        seq->seq[seq->length] = '\000';

        return (seq->ok = 1);
}

/* -------------------------------------------- */
/* lecture/redimensionnement d'une sequence au  */
/* format Fasta                                 */
/* <Lecture Indexee>                            */
/* returns : 0 -> memory error                  */
/*           1 -> read ok                       */
/* -------------------------------------------- */
int GetFastaSequence(FILE *streamin, FastaSequencePtr seq)
{
        long    readlen, buflen;
        char    *buffer, *tbuf;

        seq->ok = 0;                            /* assume error         */

        buflen = seq->length = 0; 
        
        (void) fseek(streamin, seq->offset, SEEK_SET);

        buffer = sNextIOBuffer(streamin, READ_NEXT, INDEXED);

        if (! (buffer && (*buffer == '>')))     /* sync error           */
            return 0;                           /* last sequence        */
        
        if (seq->offset)
            seq->offset -= (strlen(buffer) + 1);

        (void) strcpy(seq->name,    GetFastaName(buffer));
        
        (void) strcpy(seq->comment, GetFastaComment(buffer));
        
        while (buffer = sNextIOBuffer(streamin, READ_NEXT, INDEXED)) {

            if (*buffer == '>')
                break;

#if SECURE          
            readlen = sCountAlpha(buffer);
#else
            readlen = strlen(buffer);
#endif
        
            buflen +=  readlen;
   
            if (buflen >= seq->bufsize) {
            
                if (! (tbuf = REALLOC(char, seq->seq, 2 * buflen + 1)))
                   return 0;                 

                seq->seq = tbuf;
                
                seq->bufsize = 2 * buflen + 1;
                
            }           
#if SECURE
            sStrcpyAlpha(seq->seq + seq->length, buffer);
#else
            (void) memcpy(seq->seq + seq->length, buffer, readlen);
#endif
        
            seq->length = buflen;
        
        }

        seq->seq[seq->length] = '\000';

        return (seq->ok = 1);
}

/* -------------------------------------------- */
/* ecriture d'une sequence au format Fasta      */
/* -------------------------------------------- */
void WriteFastaSequence(FILE *streamou, FastaSequencePtr seq, 
                        int char_per_line)
{
        long  i, nlines, rest;
        char *buf, *end, tempo;

        (void) fputc('>', streamou);
        (void) fputs((*(seq->name)    ? seq->name    : "no_name"), streamou);
        (void) fputc(' ', streamou);
        (void) fputs((*(seq->comment) ? seq->comment : "\000"), streamou);
        (void) fputc(LINE_FEED, streamou);

        nlines = seq->length / char_per_line;

        buf = seq->seq;

        for (i = 0 ; i < nlines ; i++) {
            end = buf + char_per_line;
            tempo = *end;
            *end = '\000';
            (void) fputs(buf, streamou);
            (void) fputc(LINE_FEED , streamou);
            *end = tempo;
            buf += char_per_line;
        }

        if ((rest = (seq->length % char_per_line)) != 0) {
           end = buf + rest;
           tempo = *end;
           *end = '\000';
           (void) fputs(buf, streamou);
           (void) fputc(LINE_FEED , streamou);
           *end = tempo;
        }
}

/* -------------------------------------------- */
/* rewind db file                               */
/* -------------------------------------------- */
void RewindFastaDB(FILE *streamin)
{
    sRetained = 0;   /* forget previous buffer  */
    
    if (streamin)
        rewind(streamin);
}
