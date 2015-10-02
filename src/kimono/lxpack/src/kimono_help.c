/* ==================================================== */
/* @file: kim_help.c                                    */
/* Help & usage                                         */
/* @history:                                            */
/* @+   Apr. 97 <AV> first draft                        */
/* ==================================================== */

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "kimono.h"

/* ----------------------------------------------- */
/* printout help                                   */
/* ----------------------------------------------- */

#define PP fprintf(stdout, 

void PrintHelp()
{
  PP  "------------------------------------------       \n");
  PP  " Kimono Version 0.2 Nov. 05                      \n");
  PP  "------------------------------------------       \n");
  PP  "synopsis :                                       \n");
  PP  "                                                 \n");
  PP  "usage:  kimono [options] [fastafile]             \n");
  PP  "------------------------------------------       \n");
  PP  "options:                                         \n");
  PP  "----- General options   ------------------       \n");
  PP  "-a string    : use symbols from string as [A]lphabet\n");
  PP  "               default = acgt                    \n");
  PP  "                                                 \n");      
  PP  "-b string    : exclude symbols from string       \n");
  PP  "               in practice replace symbols from  \n");
  PP  "               string by X's                     \n");
  PP  "               default = <none>                  \n");
  PP  "                                                 \n");      
  PP  "-c selector  : use selector as [C]ount selector  \n");
  PP  "               selector:                         \n");
  PP  "                 count : counts symbols          \n");
  PP  "                 cumul : cumulative plot         \n");
  PP  "                 codon : count codons            \n");
  PP  "                 genes : count genes             \n");
  PP  "               default = count                   \n");
  PP  "                                                 \n");      
  PP  "-e code      : g[E]netic code                    \n");
  PP  "               code:                             \n");
  PP  "                 0 : Universal                   \n");
  PP  "                 1 : mito yeast                  \n");
  PP  "                 2 : mito vertebrate             \n");
  PP  "                 3 : filamentous fungi           \n");
  PP  "                 4 : mito insects & platyhelminthes\n");
  PP  "                 5 : Candida cylindracea         \n");
  PP  "                 6 : Ciliata                     \n");
  PP  "                 7 : Euplotes                    \n");
  PP  "                 8 : mito echinoderms            \n");
  PP  "               default = 0                       \n");
  PP  "                                                 \n");      
  PP  "-d selector  : use selector as [D]istance        \n");
  PP  "               selector:                         \n");
  PP  "                chi    : chi square distance     \n");
  PP  "                hell   : hellinger distance      \n");
  PP  "                euclid : euclidian distance      \n");
  PP  "               default = chi                     \n");
  PP  "                                                 \n");      
  PP  "-h           : [H]elp - print <this> help        \n");
  PP  "                                                 \n");      
  PP  "-i filename  : genes [I]nclude file              \n");
  PP  "               default = <none>                  \n");
  PP  "                                                 \n");      
  PP  "-m selector  : use selector as [M]ixer           \n");
  PP  "               selector:                         \n");
  PP  "                basic  : basic mixer             \n");
  PP  "                cross  : cross strands mixer     \n");
  PP  "               default = basic                   \n");
  PP  "                                                 \n");      
  PP  "-p [1[2[3]]] : [P]hase(s) for counts             \n");
  PP  "               default = <nophase>               \n");
  PP  "                                                 \n");      
  PP  "-P [1[2[3]]] : [P]hase(s) for counts             \n");
  PP  "               same as -p + base complement on   \n");
  PP  "               reverse strand			        \n");
  PP  "               default = <nophase> <nocomplement>\n");
  PP  "                                                 \n");      
  PP  "-s step      : window [S]tep in bp               \n");
  PP  "               default = 1000                    \n");
  PP  "                                                 \n");      
  PP  "-w width     : window [W]idth in bp              \n");
  PP  "               default = <see -W>                \n");
  PP  "                                                 \n");      
  PP  "-W %%width   : window [W]idth in %% of chromosome length\n");
  PP  "               default = 25%%                    \n");
  PP  "                                                 \n");      
  PP  "-x filename  : genes e[X]clude file              \n");
  PP  "               default = <none>                  \n");
  PP  "                                                 \n");      
  PP  "------------------------------------------       \n");
  PP  "file formats                                     \n");
  PP  "                                                 \n");
  PP  "fastafile: fasta format                          \n");
  PP  "include/exclude files: one gene per line :       \n");
  PP  " from to strand [comments]                       \n");
  PP  "   from, to : integers  (from <= to)             \n");
  PP  "   strand : C|D|+|-                              \n");
  PP  "--------------------------------                 \n");
  PP  " see documentation for more information          \n");
  PP  "--------------------------------                 \n");
  PP  "                                                 \n");
}

#undef PP

/* ----------------------------------------------- */
/* printout usage and exit                         */
/* ----------------------------------------------- */

#define PP fprintf(stderr, 

int ExitUsage(int stat)
{
  PP  "usage:                                           \n");
  PP  "type \"kimono -h\" for help                      \n");
  
  if (stat)
    exit(stat);
                  
  return stat;
}

#undef PP
