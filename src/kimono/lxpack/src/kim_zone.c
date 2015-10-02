/* ==================================================== */
/* @file: kim_zone.c                                    */
/* zones handling                                       */
/* @history:                                            */
/* @+   Jul. 06 <AV> first draft                        */
/* ==================================================== */

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "kimono.h"

/* -------------------------------------------- */
/* internal functions                           */
/* -------------------------------------------- */

/* -------------------------------------------- */
/* API                                          */
/* -------------------------------------------- */

/* ----------------------------------------------- */
/* free zones                                      */
/* ----------------------------------------------- */

SetOfIntervals *FreeZones(SetOfIntervals *set) {

  if (set) {
    if (set->interval)
      FREE(set->interval);
    FREE(set);
  }
  
  return NULL;
}

/* ----------------------------------------------- */
/* read zones                                      */
/* note: zones are read as 1-based intervals       */
/*       but stored as 0-based intervals           */
/* ----------------------------------------------- */

SetOfIntervals *ReadZones(int seqlength) {

  int from, to;
  SetOfIntervals *set;
  Interval *tbuf;
  char buffer[BUFSIZ], err[BUFSIZ];
  char strand[2];
  
  set = NEW(SetOfIntervals);
  set->capacity  = 256;       /* just to start with */
  set->size      = 0;
  set->interval  = NEWN(Interval, set->capacity);
  
  while (fgets(buffer, sizeof(buffer), stdin)) {
  
    StrTrimLeading(buffer);
    
    if (   (*buffer == '#')
        || (*buffer == '/')
        || (strlen(buffer) == 0))
      continue;
      
    if (sscanf(buffer, "%d%d%1s", &from, &to, strand) != 3) {
      (void) sprintf(err, "include buffer: invalid format : %s", buffer);
      Erreur(err, -1);
      // return FreeZones(set);
      continue;
    }

    if ((from > to) || (from <= 0) || (to <= 0) || (to > seqlength)) {
      (void) sprintf(err, "include buffer: invalid from/to values : %s", buffer);
      Erreur(err, -1);
      // return FreeZones(set);
      continue;
    }
      
    set->size++;
    
    if (set->size >= set->capacity) {
      set->capacity *= 2;
      if (! (tbuf = REALLOC(Interval, set->interval, set->capacity))) {
        MemoryErreur("ReadZones", 2);
        return FreeZones(set);
      }
      set->interval = tbuf;
    }
    
    tbuf = set->interval + set->size - 1;
    tbuf->from   = from - 1;  
    tbuf->to     = to - 1;
    tbuf->strand = (strcmp(strand, "D") && strcmp(strand, "+") ? -1 : +1);  
  }

  return set;  
}
