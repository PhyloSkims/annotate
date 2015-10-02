/* ==================================================== */
/* @file: kim_util.c                                    */
/* general utilities functions                          */
/* @history:                                            */
/* @+   Apr. 97 <AV> first draft                        */
/* ==================================================== */

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "kimono.h"

/* -------------------------------------------- */
/* error reporting                              */
/* -------------------------------------------- */

int Erreur(char *msg, int stat)
{

  fprintf(stderr, "*%s* [%d] %s\n", (stat > 0 ? "Error" : "Warning"), 
                                    stat, msg);

  if (stat > 0)
    exit(stat);
  
  fflush(stderr);
  
  return 0;
}

/* -------------------------------------------- */
/* memory error                                 */
/* -------------------------------------------- */

int MemoryErreur(char *where, int stat)
{
  char buffer[BUFSIZ];
  
  sprintf(buffer, "Not enough memory in %s", where);
  
  return Erreur(buffer, stat);        
}

/* -------------------------------------------- */
/* open a file and assign it to stdin           */
/* -------------------------------------------- */

int AssignToStdin(char *filename)
{
  char buffer[256];
      
  if (! freopen(filename, "r", stdin)) {
    sprintf(buffer, "cant open file %s to read", filename);
    return Erreur(buffer, -1);
  }

  return 1;
}

/* ---------------------------------------------------- */
/* Trim leading spaces from string 'str'                */
/* ---------------------------------------------------- */
char *StrTrimLeading(char *str)
{
  char *sb, *sn;
      
  if (! (str && isspace(*str)))
    return str;
              
  for (sb = sn = str ; isspace(*sn) ; sn++) 
    /* nop */ ;

  while (*sn)
    *sb++ = *sn++;

  *sb = '\000';
      
  return str;
}

