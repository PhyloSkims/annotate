/* ---------------------------------------------------------------- */
/* Copyright (c) Atelier de BioInformatique                         */
/* @file: fasta_io.c                                                */
/* @desc: sequence IO in fasta format                               */
/*                                                                  */
/* @history:                                                        */
/* @+     <Gloup> : Aug 92 : first version                          */
/* @+       <Gloup> : Nov 95 : last revised version                 */
/* ---------------------------------------------------------------- */

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "fasta_io.h"

#define SECURE   1         /* 1 : secure copy ; 0 : quickest copy  */

#define READ_NEXT 0        /* force to read next io buffer line    */
#define PUSH_BACK 1        /* force to push back io buffer line    */

#define SERIAL    1        /* serial reading                       */
#define INDEXED   0        /* indexed reading                      */

#ifdef MAC_OS_C
#define LINE_FEED '\r'
#else
#define LINE_FEED '\n'
#endif

static int sRetained = 0;  /* 1 : reuse last buffer line           */

/* -------------------------------------------- */
/* @static: lecture bufferisee                  */
/* this is the main IO function                 */
/* -------------------------------------------- */

static char * sNextIOBuffer(FILE *streamin, int retain, int serial)
{
  char   *buf, *end;
  int    reused;
  
  static char sBuffer[BUFSIZ];  /* BUFSIZ in <stdio.h>   */
  
  reused = (retain || sRetained) && serial;

  sRetained = retain;

  buf = (   reused 
          ? sBuffer 
          : fgets(sBuffer, sizeof(sBuffer), streamin));   

  if (! buf)                    /* end of file          */
    return NULL;

  end = buf + strlen(buf) - 1;

  if (*end == LINE_FEED)       /* remove trailing \n     */
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

static char *sGetFastaName(char *buffer)
{
  char c[2];
  static char sName[FASTA_NAMLEN+2];

  *c = buffer[FASTA_NAMLEN];
  buffer[FASTA_NAMLEN] = '\000';
  
  if (sscanf(buffer + 1, "%s", sName) != 1)
    (void) strcpy(sName, "<no Name>");

  buffer[FASTA_NAMLEN] = *c;
   
  return sName;
}

/* -------------------------------------------- */
/* returns sequence comment (FASTA)             */
/* -------------------------------------------- */

static char *sGetFastaComment(char *buffer)
{
  char   c[2], *space;
  static char sComment[FASTA_COMLEN+2];

  *c = buffer[FASTA_COMLEN];
  buffer[FASTA_COMLEN] = '\000';
  
  space = sNextSpace(buffer);
  
  (void) strcpy(sComment, (space ? space + 1 : "<no comment>"));

  buffer[FASTA_COMLEN] = *c;

  return sComment;
}

/* -------------------------------------------- */
/* API                                          */
/* -------------------------------------------- */

/* -------------------------------------------- */
/* liberation d'une sequence                    */
/* -------------------------------------------- */

FastaSequencePtr FreeFastaSequence(FastaSequencePtr seq)
{
  if (seq) {
    if (seq->seq)
      FREE(seq->seq);
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

/* ----------------------------------------------- */
/* lecture <Serie> d'une sequence au format  Fasta */
/* returns : 0 -> last sequence                    */
/*           1 -> more to read                     */
/*           <but> you must check seq->ok          */
/*           that may be 0 (memory error)          */
/* ----------------------------------------------- */

int ReadFastaSequence(FILE *streamin, FastaSequencePtr seq)
{
  long  readlen, buflen;
  char  *buffer, *tbuf;

  seq->ok = 0;                         /* assume error    */

  buflen = seq->length = 0; 
  
  seq->offset = ftell(streamin);

  buffer = sNextIOBuffer(streamin, READ_NEXT, SERIAL);

  if (! (buffer && (*buffer == '>')))   /* sync error     */
    return 0;                           /* last sequence  */
  
  if (seq->offset)
    seq->offset -= (strlen(buffer) + 1);

  (void) strcpy(seq->name, sGetFastaName(buffer));
  
  (void) strcpy(seq->comment, sGetFastaComment(buffer));
  
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

/* ------------------------------------------------- */
/* lecture <Indexee> d'une sequence au format  Fasta */
/* use the internal seq->offset variable to position */
/* on the correct sequence.                          */ 
/* returns : 0 -> sync or memory error               */
/*           1 -> read ok                            */
/* ------------------------------------------------- */

int GetFastaSequence(FILE *streamin, FastaSequencePtr seq)
{
  long  readlen, buflen;
  char  *buffer, *tbuf;

  seq->ok = 0;            /* assume error   */

  if (seq->offset < 0)
    return 0;                           /* bad offset    */

  buflen = seq->length = 0; 
  
  (void) fseek(streamin, seq->offset, SEEK_SET);

  buffer = sNextIOBuffer(streamin, READ_NEXT, INDEXED);

  if (! (buffer && (*buffer == '>')))   /* sync error     */
    return 0;                           /* last sequence  */
  
  (void) strcpy(seq->name, sGetFastaName(buffer));
  
  (void) strcpy(seq->comment, sGetFastaComment(buffer));
  
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
  long  i, nlines, nrest, nchar;
  char *buf, *end, tempo;

  (void) fputc('>', streamou);
  (void) fputs((*(seq->name) ? seq->name : "<no name>"), streamou);
  (void) fputc(' ', streamou);
  (void) fputs((*(seq->comment) ? seq->comment : "<no comment>"), streamou);
  (void) fputc(LINE_FEED, streamou);

  nlines = seq->length / char_per_line;
  nrest  = seq->length % char_per_line;

  if (nrest) 
    nlines++;
      
  buf = seq->seq;

  for (i = 1 ; i <= nlines ; i++) {
    nchar = ((nrest && (i == nlines)) ? nrest : char_per_line);
    end = buf + nchar;
    tempo = *end;
    *end = '\000';
    (void) fputs(buf, streamou);
    (void) fputc(LINE_FEED , streamou);
    *end = tempo;
    buf += nchar;
  }
  
}

/* -------------------------------------------- */
/* rewind db file                               */
/* -------------------------------------------- */
void RewindFastaDB(FILE *streamin)
{
  sRetained = 0; /* forget previous buffer  */
    
  if (streamin)
    rewind(streamin);
}

